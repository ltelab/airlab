clear; close all;

%% User parameters & initialisation
% User parameters
K=4; 
N_it=1; 
method='svm';
type_classif='multiclass'; 
use_cost_weights = false;
target = 'svm';
normalization_type = 'standardization';
dynamic_feat_transfo = false;
random_CV = false;

% Path to training data
if strcmp(target,'2DS') || strcmp(target,'svm')
    dir_data = '../training_set/2DS_smooth';
elseif strcmp(target,'HVPS')
    dir_data = '../training_set/HVPS';
end
data_filenames = dir(fullfile(dir_data,'*.mat'));
data_filenames = {data_filenames.name}';
data_picnames = dir(fullfile(dir_data,'*.png'));
data_picnames = {data_picnames.name}';

% Time interval 
t_str_start = '20150101000000';
t_str_stop  = '20180101000000';

% chose feat_vec here
load('features_opt_4fold_10it_alpha0.0001_lambda0.01_2DS_3500samples.mat');
n_desc = 17;
feat_vec = feat_mat(:,n_desc+1);
feat_vec(feat_vec==0) = [];
feat_vec = [feat_vec];
% add ratio touching the boarder ?
icpca_feats = [70:79]';
feat_vec = [feat_vec; icpca_feats];
%feat_vec(end+1) = 68;
%feat_vec(end+1) = 69;

% hyperparameters for ML method
if strcmp(target,'2DS')
    parameters_method =  {0.0001,0.01,10000,0,10000}; % 0.0001 or 0.001 for stepsize / 1 or 0.1 for lambda
elseif strcmp(target,'HVPS')
    parameters_method = {0.001,1,10000,0,10000};
elseif strcmp(target,'riming')
    parameters_method = {0.0001,0.01,1000,0,10000};
elseif strcmp(target,'melting')
    parameters_method = {0.001,0.01,1000,0,10000};
elseif strcmp(target,'svm')
    parameters_method = {100,0.01,'rbf'}; % 100 / 0.01 for 17 desc || 1 0.01 for all
end



%% Data loading

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
    yR(idx_unknown) = [];
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

%% Classification 
N = size(X,1);
D = size(X,2);
N_classes = length(unique(y));
feat_vec_dummy = 1:1:size(X,2);
Cout = CrossValidation(method,type_classif,parameters_method,K,N_it,X,y,random_CV,1,1,use_cost_weights,dynamic_feat_transfo,feat_vec_dummy);

    