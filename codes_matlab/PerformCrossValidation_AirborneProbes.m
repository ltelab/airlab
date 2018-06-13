clearvars; close all;

%% User parameters & initialisation
% User parameters
K=4; 
N_it=10; 
method='logistic'; 
use_cost_weights = true;
target_type = 'subclassif'; % main | subclassif
target = 'AG-BR_AG-O';
normalization_type = 'standardization';
dynamic_feat_transfo = false;
random_CV = false;
grid_search = false;
disp_summary = true;

if strcmp(target_type,'main')
    type_classif = 'multiclass'; 
else
    type_classif = 'binary';
end


% Path to training data
if strcmp(target,'2DS') || strcmp(target,'svm')
    dir_data = '../training_set/2DS_rospla2_final';
elseif strcmp(target,'HVPS')
    dir_data = '../training_set/HVPS_new3';
elseif strcmp(target,'CPI')
    dir_data = '../training_set/CPI_smooth0_icpca0';
elseif strcmp(target,'AG-BR_AG-O')
    dir_data = '../training_set/subclassification/bulrosagg_vs_others/all';
end
data_filenames = dir(fullfile(dir_data,'*.mat'));
data_filenames = {data_filenames.name}';
%data_picnames = dir(fullfile(dir_data,'*.png'));
%data_picnames = {data_picnames.name}';

% Time interval 
t_str_start = '20150101000000';
t_str_stop  = '20180101000000';

% chose feat_vec here
if strcmp(target_type,'main')
    if strcmp(target,'2DS')
        load('feat_opt/BER_based/2DS_4217s_10it_4k.mat');
    elseif strcmp(target,'HVPS')
        load('feat_opt/BER_based/HVPS_1426s_10it_4k.mat');
    elseif strcmp(target,'CPI')
        load('feat_opt/BER_based/CPI_2964s_10it_4k.mat');
    else
        fprintf('Error : target %s not reckognized ! \n',target);
    end
    
    n_desc = 15;
    feat_vec = feat_mat(:,n_desc+1);
    feat_vec(feat_vec==0) = [];
      
elseif strcmp(target_type,'subclassif')

    n_desc = 2;
    if strcmp(target,'AG-BR_AG-O')
        feat_vec = [83 33];
    else
        fprintf('Error : target %s not reckognized ! \n',target);
    end
    
end


%feat_vec = [feat_vec; [99:111]'];
% add ratio touching the boarder ?
%icpca_feats = []';
%feat_vec = [feat_vec; icpca_feats];
%feat_vec(end+1) = 68;
%feat_vec(end+1) = 69;

% hyperparameters for ML method
if strcmp(method,'svm')
    parameters_method = {100,0.01,'rbf'}; % 100 / 0.01 for 17 desc || 1 0.01 for all
elseif strcmp(target,'2DS')
    parameters_method =  {0.0001,0.1,5000,0,5000}; % 0.0001 or 0.001 for stepsize / 1 or 0.1 for lambda // after grid search with 20 features : {0.00016681,0.22,5000,0,5000} // fine tuning : {0.00014384,0.38,5000,0,5000}
elseif strcmp(target,'HVPS')
    parameters_method = {0.001,1,5000,0,5000};
elseif strcmp(target,'CPI')
    parameters_method = {0.0001,0.1,5000,0,5000}; % old educated guess {0.0001,0.1,5000,0,5000}; // fine tuning : {4.8329e-4,0.1833,5000,0,5000}
elseif strcmp(target,'riming')
    parameters_method = {0.0001,0.01,1000,0,10000};
elseif strcmp(target,'melting')
    parameters_method = {0.001,0.01,1000,0,10000};
elseif strcmp(target,'AG-BR_AG-O')  
    parameters_method = {0.0001,0.1,5000,0,5000};
end



%% Data loading

% Load the training matrix X
[X,Xlab,Xname,Xt] = load_processed_2DS_data(dir_data,t_str_start,t_str_stop,feat_vec);
X_ini = X;

% Load the labels vector y
if strcmp(target_type,'main')
    y = load_2DS_labels(dir_data,t_str_start,t_str_stop);
elseif strcmp(target_type,'subclassif')
    [~,y] = load_2DS_labels(dir_data,t_str_start,t_str_stop,1);
end

% remove NaN-values in X,y
idx = find(isnan(sum(X,2)));
X(idx,:) = [];
y(idx) = [];

% remove truncated particles
% idx = find(X(:,15) > 0.325 & y < 5);
% X(idx,:) = [];
% y(idx) = [];

% remove unwanted labels
idx_unknown = find(y<0 | y>7);
if ~isempty(idx_unknown)
    y(idx_unknown) = [];
    X(idx_unknown,:) = [];
    %data_picnames(idx_unknown) = [];
    data_filenames(idx_unknown) = [];
    fprintf('%u samples discarded because not labelled correctly \n',length(idx_unknown));
end

% load categories names
if strcmp(target,'2DS')
    labels = {'AG','CC','CP','BR','QS','OT','PC'};
elseif strcmp(target,'HVPS')
    labels = {'AG','CC','CP','BR','QS'};
elseif strcmp(target,'CPI')
    labels = {'AG','CC','CP','BR','QS','PC'};
elseif strcmp(target,'AG-BR_AG-O')
    labels = {'AG-BR','AG-O'};
end

% summary
if disp_summary
   fprintf('\n\n*****************************\n');
   for i=1:numel(unique(y))
       
       fprintf('Class %u (%s) : %5u samples (%2.1f%%) \n',i,labels{i},sum(y==i),100*sum(y==i)/numel(y));
       
   end
   fprintf('*****************************\n'); 
end


%% Features transformation
if ~dynamic_feat_transfo

D = size(X,2);

for i=1:D
      
    skew = skewness(X(:,i));
    % save the skewness in model
    classif_params.skew(i,1) = skew;
    
    if skew > 1
        fprintf('Feat %u, skewness = %2.2f \n',i,skew)
        X(:,i) = log(abs(X(:,i)+1)); 
        
    elseif skew > 0.75
        X(:,i) = sqrt(abs(X(:,i)));
        
    elseif skew < -1
        X(:,i) = exp(X(:,i));
        
    elseif skew < -0.75
        X(:,i) = X(:,i).^2;
        
    end

    if strcmp(normalization_type,'standardization') 
        
        tmp_mean = mean(X(:,i));
        tmp_std  = std(X(:,i));
        if tmp_std == 0
            tmp_std = 1;
        end
        
        if i<1000000000
            X(:,i) = (X(:,i) - tmp_mean)/tmp_std;
        else
            X(:,i) = (X(:,i) - tmp_mean)/(5*tmp_std);
        end   
        
        classif_params.mean(i,1) = tmp_mean;
        classif_params.std(i,1) = tmp_std;
        
    elseif strcmp(normalization_type,'rescaling')
        
        tmp_min = min(X(:,i));
        tmp_max = max(X(:,i));
        X(:,i) = (X(:,i) - tmp_min)/(tmp_max - tmp_min);
        
    end
    
end

end

%% visual illustration of norm + Gauss
if 1
    disp_feat = 2;
    figure; subplot(121); hold on; box on;
    histogram(X_ini(:,disp_feat),20,'Normalization','count','facecolor',[0.3 0.3 0.3]);
    xlabel('Before');
    ylabel('Count');
    title('Feature #14');
    subplot(122); hold on; box on;
    histogram(X(:,disp_feat),20,'Normalization','count','facecolor',[0.3 0.3 0.3]);
    xlabel('After');
    
end


%% Classification 
if ~grid_search
    
    N = size(X,1);
    D = size(X,2);
    N_classes = length(unique(y));
    feat_vec_dummy = 1:1:size(X,2);
    Cout = CrossValidation(method,type_classif,parameters_method,K,N_it,X,y,random_CV,1,1,use_cost_weights,dynamic_feat_transfo,feat_vec,labels);
    

%% hyperparameters space grid search
else

    %{0.00016681,0.22,5000,0,5000}
    
    N = size(X,1);
    D = size(X,2);
    N_classes = length(unique(y));
    feat_vec_dummy = 1:1:size(X,2);

    step_vec = logspace(-5,-3,20);
    lambda_vec = logspace(-2,1,20);
    alpha_vec = [0];
    it_vec = [5000 10000];
    
    BER_Te = zeros(numel(step_vec),numel(lambda_vec));
    BER_Tr = zeros(numel(step_vec),numel(lambda_vec));
    kappa_Te = zeros(numel(step_vec),numel(lambda_vec));
    kappa_Tr = zeros(numel(step_vec),numel(lambda_vec));
    OA_Te = zeros(numel(step_vec),numel(lambda_vec));
    OA_Tr = zeros(numel(step_vec),numel(lambda_vec));
    
    count = 0;
    count_tot = numel(step_vec)*numel(lambda_vec)*numel(alpha_vec)*numel(it_vec);
    
    for i=1:numel(step_vec)
        parfor j=1:numel(lambda_vec)
            %for k=1:numel(it_vec)
                %for l=1:numel(alpha_vec)
                    
                    %count = count + 1;
                    %fprintf('%u / %u \n',count,count_tot);
                    fprintf('%u/%u || %u/%u || %u/%u || %u/%u \n',i,numel(step_vec),j,numel(lambda_vec),1,numel(it_vec),1,numel(alpha_vec));
                    
                    parameters_method =  {step_vec(i),lambda_vec(j),5000,alpha_vec(1),it_vec(1)};
                    Cout = CrossValidation(method,type_classif,parameters_method,K,N_it,X,y,random_CV,0,0,use_cost_weights,dynamic_feat_transfo,feat_vec_dummy);
                    BER_Te(i,j) = Cout.BER_Te;
                    BER_Tr(i,j) = Cout.BER_Tr;
                    kappa_Te(i,j) = Cout.kappa_Te;
                    kappa_Tr(i,j) = Cout.kappa_Tr;
                    OA_Te(i,j) = Cout.OA_Te;
                    OA_Tr(i,j) = Cout.OA_Tr;  
                                        
                %end
            %end
        end
    end
    
end
    
%%
if 0

kappa_max = max(kappa_Te(:));
fprintf('%2.2f%% \n',100*kappa_max);
n = 1;

for i=1:numel(step_vec)
    for j=1:numel(lambda_vec)
        for k=1:numel(it_vec)
            for l=1:numel(alpha_vec)
                
                if kappa_Te(i,j,l,k) == kappa_max
                    
                    idx(n,:) == [i j l k];
                    fprintf('%5.8f %5.2f %2.2f %u \n',step_vec(i),lambda_vec(j),alpha_vec(l),it_vec(k));
                    n = n+1;
                    
                end
                
            end
        end
    end
end

end


    


    