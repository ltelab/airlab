clearvars; close all;


K=10; 
N_it=20; 
method='logistic';
type_classif='multiclass'; 
use_cost_weights = true;
target = '2DS';
normalization_type = 'standardization';
dynamic_feat_transfo = true; %<- this option is useless here as dynamic feature transformation is always enabled

% useful parameters?
features_ranking = false;
display_confmat = true;
label_version = '1.1';


if strcmp(target,'2DS')
    parameters_method =  {0.00014384,0.1,5000,0,5000}; % 0.0001 or 0.001 for stepsize / 1 or 0.1 for lambda
elseif strcmp(target,'CPI')
    parameters_method = {0.0001,0,5000,0,5000};
elseif strcmp(target,'HVPS')
    parameters_method = {0.001,1,1000,0,10000};
elseif strcmp(target,'riming')
    parameters_method = {0.0001,0.01,1000,0,10000};
elseif strcmp(target,'melting')
    parameters_method = {0.001,0.01,1000,0,10000};
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% path to training data
if strcmp(target,'2DS')
    dir_data = '../training_set/2DS_4000_smooth0_icpca1/mat';
elseif strcmp(target,'HVPS')
    dir_data = '/home/praz/Documents/Cloud_Probes/training_set/HVPS';
elseif strcmp(target,'CPI')
    dir_data = '/home/praz/Documents/airlab/training_set/CPI12';
end
data_filenames = dir(fullfile(dir_data,'*.mat'));
data_filenames = {data_filenames.name}';
data_picnames = dir(fullfile(dir_data,'*.png'));
data_picnames = {data_picnames.name}';

% Time interval 
t_str_start = '20150101000000';
t_str_stop  = '20180101000000';

% chose feat_vec here
if strcmp(target,'2DS')
    load('feat_opt/2DS/rand_4fold_10it_3990N_98D.mat');
elseif strcmp(target,'HVPS')
    load('feat_opt/features_opt_4fold_20it_alpha0.001_lambda1_HVPS.mat');
elseif strcmp(target,'CPI')
    load('feat_opt/CPI/rand_4fold_10it_more_samples.mat');
else
    fprintf('Error : target %s not reckognized ! \n',target);
end
n_desc = 15;
feat_vec = feat_mat(:,n_desc+1);
feat_vec(feat_vec==0) = [];
feat_vec = [feat_vec];

% add ratio touching the boarder ? 
%feat_vec(end+1) = 67;
%feat_vec(end+1) = 68;
%feat_vec(end+1) = 69;

% Load the training matrix X
[X,Xlab,Xname,Xt] = load_processed_2DS_data(dir_data,t_str_start,t_str_stop,feat_vec);

% Load the labels vector y
y = load_2DS_labels(dir_data,t_str_start,t_str_stop);

% remove NaN-values in X,y
idx = find(isnan(sum(X,2)));
X(idx,:) = [];
y(idx) = [];

% remove unknown labels
idx_unknown = find(y<=0);
if ~isempty(idx_unknown)
    y(idx_unknown) = [];
    X(idx_unknown,:) = [];
    data_picnames(idx_unknown) = [];
    data_filenames(idx_unknown) = [];
    fprintf('%u samples discarded because not labelled correctly \n',length(idx_unknown));
end

% load categories names
if strcmp(target,'2DS')
    labels = {'Agg','Col','Gra','Ros','Sph','Oth'};
elseif strcmp(target,'HVPS')
    labels = {'Agg','Col','Gra','Ros','Sph'};
elseif strcmp(target,'CPI')
    labels = {'Agg','Col','Gra','Ros','Sph','Pla'};
end

y = real(y);


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

D = size(X,2);

% for i=1:D
%       
%     skew = skewness(X(:,i));
%     % save the skewness in model
%     classif_params.skew(i,1) = skew;
%     
%     if skew > 1
%         
%         X(:,i) = log(abs(X(:,i)+1));
%         
%     elseif skew > 0.75
%         
%         X(:,i) = sqrt(abs(X(:,i)));
%         
%     elseif skew < -1
%         
%         X(:,i) = exp(X(:,i));
%         
%     elseif skew < -0.75
%         
%         X(:,i) = X(:,i).^2;
%         
%     end
% 
%     if strcmp(normalization_type,'standardization') 
%         
%         tmp_mean = mean(X(:,i));
%         tmp_std  = std(X(:,i));
%         if i<1000000000
%             X(:,i) = (X(:,i) - tmp_mean)/tmp_std;
%         else
%             X(:,i) = (X(:,i) - tmp_mean)/(5*tmp_std);
%         end   
%         
%         classif_params.mean(i,1) = tmp_mean;
%         classif_params.std(i,1) = tmp_std;
%         
%     elseif strcmp(normalization_type,'rescaling')
%         
%         tmp_min = min(X(:,i));
%         tmp_max = max(X(:,i));
%         X(:,i) = (X(:,i) - tmp_min)/(tmp_max - tmp_min);
%         
%     end
%     
% end


%% Algorithm training & testing

% split between train and test
setSeed(100);
idx_perm = randperm(N);
Nk = floor(N/K);
idxTe = idx_perm(1:Nk);
idxTr = idx_perm(Nk+1:end);
yTe = y(idxTe);
XTe_ini = X(idxTe,:);
yTr_tot = y(idxTr);
XTr_tot = X(idxTr,:);
%train_sampling = [0.05 0.10 0.15 0.20 0.25 0.3 0.35 0.4 0.45 0.5 0.55 0.6 0.65 0.7 0.75 0.8 0.85 0.9 0.95];
train_sampling = [0.02:0.02:1];
%[0.02 0.03 0.04 0.05 0.06 0.07 0.08 0.09 0.10 0.15 0.20 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1];

for i=1:N_it
    
    fprintf('\n Iteration %u/%u : \n\n',i,N_it);

    idx_perm = randperm(length(yTr_tot));
    for k=1:length(train_sampling)

        idx2keep =idx_perm(1:floor(train_sampling(k)*length(yTr_tot)));
        yTr = yTr_tot(idx2keep);
        XTr = XTr_tot(idx2keep,:);
        
        XTe = XTe_ini;
        
        % Normalization based on training set only
        for d=1:D

            skew = skewness(XTr(:,d));

            if skew > 1

                XTr(:,d) = log(abs(XTr(:,d)+1));
                XTe(:,d) = log(abs(XTe(:,d)+1));

            elseif skew > 0.75

                XTr(:,d) = sqrt(abs(XTr(:,d)));
                XTe(:,d) = sqrt(abs(XTe(:,d)));

            elseif skew < -1

                XTr(:,d) = exp(XTr(:,d));
                XTe(:,d) = exp(XTe(:,d));

            elseif skew < -0.75

                XTr(:,d) = XTr(:,d).^2;
                XTe(:,d) = XTe(:,d).^2;

            end

            if strcmp(normalization_type,'standardization') 

                tmp_mean = mean(XTr(:,d));
                tmp_std  = std(XTr(:,d));
                
                XTr(:,d) = (XTr(:,d) - tmp_mean)/tmp_std;
                XTe(:,d) = (XTe(:,d) - tmp_mean)/tmp_std;
                
            end
 
        end

        % Fit model using the wrapping function
        model = trainClassifier(method, type_classif, yTr, XTr, parameters_method);
        
%         % Features ranking
%         betas = abs(model);
%         beta_sum = sum(betas,2);
%         beta_mean = mean(betas,2);
%         [~,beta_rank_sum] = sort(beta_sum,'descend');
%         [~,beta_rank_mean] = sort(beta_mean,'descend');
%         beta_rankvec_sum = [beta_rankvec_sum beta_rank_sum];
%         beta_rankvec_mean = [beta_rankvec_mean beta_rank_mean];
        
        
        % Predict
        [predTe,scoresTe] = predictClassifier(method,type_classif,model,XTe);
        [predTr,scoresTr] = predictClassifier(method,type_classif,model,XTr);  

        % Compute accuracy
         BER_Te(i,k) = computeBER(yTe,predTe);
         BER_Tr(i,k) = computeBER(yTr,predTr);
         OA_Te(i,k) = computeOA(yTe,predTe);
         OA_Tr(i,k) = computeOA(yTr,predTr);
         kappa_Te(i,k) = computeKAPPA(yTe,predTe);
         kappa_Tr(i,k) = computeKAPPA(yTr,predTr);

         fprintf('size of training sample : %2.2f%%',train_sampling(k)*100);
         fprintf('Testing BER: %.2f%%\n', BER_Te(i,k) * 100 );
         fprintf('Training BER: %.2f%%\n\n', BER_Tr(i,k) * 100 );

    end
    
end

        
%% diagnostic
figure;
plot_learning_curves(train_sampling',kappa_Te','Test',kappa_Tr','Train','% of training sample used','Kappa');
set(gca,'Fontsize',14);
figure;
plot_learning_curves(train_sampling',BER_Te','Test',BER_Tr','Train','% of training sample used','BER');
set(gca,'Fontsize',14);

