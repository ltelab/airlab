clearvars ; close all;

                % alpha = 0.0001, lambda = 0.01, mini-batch size = 100
                % momentum = 0.9, Max no of iterations = 5000
                %parameters={0.00001,0.01,100,0.9,5000}; 

% USER DEFINED PARAMETERS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

method='logistic';
type_classif='multiclass'; 
use_cost_weights = true;
target = '2DS';
normalization_type = 'standardization';
dynamic_feat_transfo = false;
parameters_method = {0.0001,0.01,10000,0,5000};

features_ranking = true;
display_confmat = true;
label_version = '1.1';


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% path to training data
dir_data = '/home/praz/Documents/airlab/training_set/2DS';

data_filenames = dir(fullfile(dir_data,'*.mat'));
data_filenames = {data_filenames.name}';
data_picnames = dir(fullfile(dir_data,'*.png'));
data_picnames = {data_picnames.name}';

% time interval to take into consideration
t_str_start = '20150101000000';
t_str_stop  = '20180101000000';

% chose feat_vec here
load('features_opt_4fold_10it_alpha0.0001_lambda0.01_2DS_3500samples.mat');
n_desc = 17;
feat_vec = feat_mat(:,n_desc+1);
feat_vec(feat_vec==0) = [];
feat_vec = [feat_vec];
%feat_vec(end+1) = 67;

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
    fprintf('%u samples discarded because labelled as "undetermined" or "ambiguous" \n',length(idx_unknown));
end

% load categories names
if strcmp(target,'2DS')
    labels = {'Agg','Col','Gra','Ros','Sph','Oth'};
elseif strcmp(target,'HVPS')
    labels = {'Agg','Col','Gra','Ros','Sph'};
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



