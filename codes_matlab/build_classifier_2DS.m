clearvars ; close all;

                % alpha = 0.0001, lambda = 0.01, mini-batch size = 100
                % momentum = 0.9, Max no of iterations = 5000
                %parameters={0.00001,0.01,100,0.9,5000}; 

% USER DEFINED PARAMETERS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

method='logistic';
use_cost_weights = true;
target_type = 'subclassif'; 
target = 'AG-BR_AG-O';
normalization_type = 'standardization';
dynamic_feat_transfo = false;

features_ranking = true;
display_confmat = true;
label_version = '1.1';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if strcmp(target_type,'main')
    type_classif = 'multiclass'; 
else
    type_classif = 'binary';
end

% path to training data
dir_data = '../training_set/subclassification/bulrosagg_vs_others/all';

data_filenames = dir(fullfile(dir_data,'*.mat'));
data_filenames = {data_filenames.name}';
data_picnames = dir(fullfile(dir_data,'*.png'));
data_picnames = {data_picnames.name}';

% time interval to take into consideration
t_str_start = '20150101000000';
t_str_stop  = '20180101000000';

% choose parameters_method
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

%% Data loading

% Load the training matrix X
[X,Xlab,Xname,Xt] = load_processed_2DS_data(dir_data,t_str_start,t_str_stop,feat_vec);

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

% remove unknown labels
idx_unknown = find(y<0);
if ~isempty(idx_unknown)
    y(idx_unknown) = [];
    X(idx_unknown,:) = [];
    data_picnames(idx_unknown) = [];
    data_filenames(idx_unknown) = [];
    fprintf('%u samples discarded because labelled as "undetermined" or "ambiguous" \n',length(idx_unknown));
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

N = size(X,1);
D = size(X,2);
N_classes = length(unique(y));

% computing the weights, if we want to penalize the cost
if use_cost_weights
    if strcmp(type_classif,'binary')
        for i=1:N_classes
            Nc(i,1) = sum(y==(i-1));
        end
    else
        for i=1:N_classes
            Nc(i,1) = sum(y==i);
        end
    end
    [Ncmax,~] = max(Nc);
    for i=1:N_classes
        weights(i,1) = Ncmax/Nc(i);
    end
    parameters_method{end+1} = weights;

end



%% Features transformation
if ~dynamic_feat_transfo

D = size(X,2);

for i=1:D
      
    skew = skewness(X(:,i));
    % save the skewness in model
    classif_params.skew(i,1) = skew;
    
    if skew > 1
        
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
  
% Initialize global confusion matrix
global_confmat_Te = zeros(N_classes);
global_confmat_Tr = zeros(N_classes);

% Initialize beta ranking vectors
beta_sum_mat = [];

% Building classifier: XTr = Xtot
XTr = X;
yTr = y;

% Fit model using the wrapping function
model = trainClassifier(method, type_classif, yTr, XTr, parameters_method);
        
% Features ranking
if features_ranking
    betas = abs(model);
    beta_sum = sum(betas,2);
    beta_sum = beta_sum(2:end,:);
    beta_sum_mat(:,end+1) = beta_sum;
end
               
% Predict
[predTr,scoresTr] = predictClassifier(method,type_classif,model,XTr);  
        
% Compute accuracy
BER_Tr = computeBER(yTr,predTr);
OA_Tr = computeOA(yTr,predTr);
kappa_Tr = computeKAPPA(yTr,predTr);

fprintf('Training BER: %.2f%%\n\n', BER_Tr * 100 );
fprintf('Training OA: %.2f%%\n\n', OA_Tr * 100 );
fprintf('Training kappa: %.2f%%\n\n', kappa_Tr * 100 );

% Confusion matrix
if display_confmat
    global_confmat = confusionmat(yTr,predTr);   
    global_confmat = round(global_confmat./sum(global_confmat(:)) * 100,2);
    figure;
    heatmap(global_confmat,labels,labels,1,'Colormap','red','ShowAllTicks',1,'UseLogColorMap',true,'Colorbar',true);
    title('Train confusion matrix');
end



classifier.type = method;
classifier.type_classif = type_classif;
classifier.parameters_method = parameters_method;
classifier.normalization = 'standardization';
classifier.model = model;
classifier.N = N;
classifier.N_classes = N_classes;
classifier.N_labels = labels;
classifier.D = D;
classifier.Xlab = Xlab;
classifier.normalization_params = classif_params;
classifier.BER = BER_Tr;
classifier.OA = OA_Tr;
classifier.kappa = kappa_Tr;
classifier.feat_vec = feat_vec;



