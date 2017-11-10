% function to predict snow habit for a series of images contained in a folder
% exemple:
% dir_data = '/home/praz/Documents/Cloud_Probes_test_set/2DS';
% t_str_start = '20150101000000';
% t_str_stop = '20170101000000';
% classifier = 'logit_test.mat';
% the blurry threshold was used for MASC images. In the case of airborne
% probes images (2DS, HVPS, CPI), just ignore it
% truncated threshold for 2DS (CPI?) = 0.325 // CPI frame_fraction = 0.3

function predict_snow_habit(dir_data,t_str_start,t_str_stop,classifier,save_results,save_validation_files,truncated_threshold,blurry_threshold)
 
    tic;

    if nargin == 6
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
            % check roi.status
            if strcmp(roi.status,'noisy')
                roi.label_ID = -2;
                roi.label_name = 'Noisy';
                
            elseif strcmp(roi.status,'too small')
                roi.label_ID = 7;
                roi.label_name = 'Small';
            
            % check roi.truncated_threshold
            elseif roi.frame_fraction >= truncated_threshold && pred(i) ~= 5 % pred(i)==5 is for quasi-spheres, which are fairly recognized even if truncated
                roi.label_ID = 8;
                roi.label_name = 'Trunc';
   
            % check classif output validity
            elseif isnan(scores(i,1))
                roi.label_ID = -1;
                roi.label_name = 'Error';
                
            else
                roi.label_ID = pred(i);
                roi.label_name = prediction_scheme{roi.label_ID};
            end
            
            roi.label_probs = scores(i,:);
            N_tot = N_tot + 1;
            save(fullfile(dir_data,Xname{i}),'roi');
             
            if save_validation_files
                if ~exist(fullfile(dir_data,'valid'),'dir')
                    mkdir(fullfile(dir_data,'valid'));
                end
                fig = figure('Visible','off'); 
                subplot(1,2,1); imshow(roi.raw_data); title(roi.label_name); set(gca,'fontsize',9);
                [p_max,idx_max] = max(roi.label_probs);
                subplot(1,2,2); imshow(roi.data_gs); title(sprintf('ID_{max} : %u, p_{max}(%u) : %2.1f',idx_max,idx_max,p_max)); set(gca,'fontsize',9);
                figname = roi.name;
                saveas(fig,fullfile(dir_data,'valid',figname),'png');   
            end
            
            clear roi;
            close all hidden;
            
        end

        total_time = toc;
        fprintf('   Done!\n');
        fprintf('***** %u images classified and saved in %2.2f seconds*****\n\n\n',N_tot,total_time);

    else 

        toc;

    end

end


