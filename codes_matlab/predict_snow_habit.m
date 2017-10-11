% function to predict snow habit for a series of images contained in a folder
% exemple:
% dir_data = '/home/praz/Documents/Cloud_Probes_test_set/2DS';
% t_str_start = '20150101000000';
% t_str_stop = '20170101000000';
% classifier = 'logit_test.mat';
% the blurry threshold was used for MASC images. In the case of airborne
% probes images (2DS, HVPS), just ignore it

function predict_snow_habit(dir_data,t_str_start,t_str_stop,classifier,save_results,blurry_threshold)
 
    tic;

    if nargin == 5
        blurry_threshold = 0;
    end

    load(classifier);
    normalization_type = classifier.normalization;
    prediction_scheme = classifier.N_labels;

    fprintf('Classifying snowflakes type in : %s...\n',dir_data);
    %% load flakes descriptors

    % load the training matrix X
    [X,Xlab,Xname,~] = load_processed_2DS_data(dir_data,t_str_start,t_str_stop,classifier.feat_vec);
    
    % if X empty, we stop
    if isempty(X)
        fprintf('Error: unable to load the matrix X of descriptors \n');
        return;
    end

    %% Features transformation
    % according to classif_params

    fprintf('Assembling processed data in X matrix and transforming features...');

    D = size(X,2);

    for i=1:D

        if classifier.normalization_params.skew(i) > 1

            X(:,i) = log(abs(X(:,i)+1));

        elseif classifier.normalization_params.skew(i) > 0.75

            X(:,i) = sqrt(abs(X(:,i)));

        elseif classifier.normalization_params.skew(i) < -1

            X(:,i) = exp(X(:,i));

        elseif classifier.normalization_params.skew(i) < -0.75

            X(:,i) = X(:,i).^2;

        end

        if strcmp(normalization_type,'standardization') 

            X(:,i) = (X(:,i) - classifier.normalization_params.mean(i))/classifier.normalization_params.std(i);

        elseif strcmp(normalization_type,'rescaling')

            disp('warning: rescaling not implemented yet!');

        end

    end

    fprintf('  Done!\n'); 

    %% Prediction
    fprintf('Classifying images according to %s classifier...',classifier.type);  
    [pred,scores] = predictClassifier(classifier.type,classifier.type_classif,classifier.model,X);    
    fprintf('   Done!\n');

    

    %% saving
    if save_results

        fprintf('Saving results in roi structures...');
        N_tot = 0;

        for i=1:numel(Xname)
            load(fullfile(dir_data,Xname{i}));
            % if the classification failed (probably because the algorithm could not retrieve some descriptors), we classify the particle as "small/other"
            if isnan(scores(i,1))
                roi.label_ID = 6;
            else
                roi.label_ID = pred(i);
            end
            roi.label_name = prediction_scheme{roi.label_ID};
            roi.label_probs = scores(i,:);
            N_tot = N_tot + 1;
            save(fullfile(dir_data,Xname{i}),'roi');
        end

        total_time = toc;
        fprintf('   Done!\n');
        fprintf('***** %u images classified and saved in %2.2f seconds*****\n\n\n',N_tot,total_time);

    else 

        toc;

    end

end


