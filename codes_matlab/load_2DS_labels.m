function [y,y2] = load_2DS_labels(dirname,t_str_start,t_str_stop,load_subclassif)

    fprintf('Load associated labels...');
    
    if exist('t_str_start','var')
        tmin = datenum(t_str_start,'yyyymmddHHMMSS');
    else
        tmin = datenum('20000101000000','yyyymmddHHMMSS');
    end
    
    if exist('t_str_stop','var')
        tmax = datenum(t_str_stop,'yyyymmddHHMMSS');
    else
        tmax = datenum('21000101000000','yyyymmddHHMMSS');
    end
    
    if ~exist('load_subclassif','var')
        load_subclassif = false;
    end

    file_list = dir(fullfile(dirname,'*.mat'));
    file_only_list = {file_list.name};
    file_list = fullfile(dirname,file_only_list);

    % loop over the snowflakes and save each feature in a vector of the main
    % structure data
    n_flakes = 0;
    y = [];
    if load_subclassif
        y2 = [];
    else
        y2 = NaN;
    end
    i = 1;
    for k=1:length(file_list)

        load(file_list{k});
        % load only pictures within the time interval
        if (~isfield(roi,'tnum') || (roi.tnum >= tmin && roi.tnum <= tmax)) %&& ~roi.is2small

            n_flakes = n_flakes + 1;
            y(i,1) = roi.label_ID; 
            if load_subclassif
                y2(i,1) = roi.label_ID2;
            end
            i = i+1;

        end

    end
    
    fprintf('   Done! %u labels loaded. \n',n_flakes);
    
end