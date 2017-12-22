% function to predict snow habit for a series of images contained in a folder
% exemple:
% dir_data = 'test_sample/2DS/Proc_images/mat';
% probe = '2DS'
% save_results = true
% save_validation_files = true
% truncated_threshold = 0.3
function predict_snow_habit(dir_data,probe,save_results,save_validation_files,truncated_threshold)
 
    tic;
    
    if nargin < 5
        truncated_threshold = 0.3;
    end

    if strcmp(probe,'2DS')
        classifier = 'logit_trained_models/logit_2DS_6classes_20171124.mat';
    elseif strcmp(probe,'HVPS')
        classifier = 'logit_trained_models/logit_HVPS_5classes_20170619.mat';
        disp('Warning : HVPS classifier was trained on a very small dataset. Classification is uncertain.');
    elseif strcmp(probe,'CPI')
        classifier = 'logit_trained_models/logit_CPI_6classes_20171124.mat';       
    end
    
    load(classifier);
    normalization_type = classifier.normalization;
    prediction_scheme = classifier.N_labels;

    fprintf('Classifying snowflakes type in : %s...\n',dir_data);
    %% load flakes descriptors

    % load the training matrix X
    [X,Xlab,Xname,~] = load_processed_2DS_data(dir_data,'','',classifier.feat_vec);
    
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
                if isfield(roi,'data_gs') % CPI
                    subplot(1,2,1); imshow(roi.raw_data); 
                    if ~isfield(roi,'noise_thresh_ini')
                        roi.noise_thresh_ini = 0;
                        roi.noise_thresh_eff = 0;
                        roi.mean_intens = 0;
                    end  
                    title(sprintf('%s, n_i=%1.2f, n_f=%1.2f, B=%1.2f',roi.label_name,roi.noise_thresh_ini,roi.noise_thresh_eff,roi.mean_intens)); set(gca,'fontsize',9);
                    [p_max,idx_max] = max(roi.label_probs);
                    if ~isfield(roi,'frame_fraction')
                        roi.frame_fraction = 0;
                    end
                    subplot(1,2,2); imshow(roi.data_gs); title(sprintf('ID_{max}:%u, p_{max}(%u):%2.1f, frame frac=%1.2f',idx_max,idx_max,p_max,roi.frame_fraction)); set(gca,'fontsize',9);
  
                else
                    [p_max,~] = max(roi.label_probs);
                    imshow(~roi.raw_data); title(sprintf('habit : %s -- prob = %2.1f',roi.label_name,p_max)); set(gca,'fontsize',9);
                    %subplot(1,2,2); imshow(~roi.data); title(sprintf('ID_{max} : %u, p_{max}(%u) : %2.1f',idx_max,idx_max,p_max)); set(gca,'fontsize',9);
                    
                end
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


