function [X,Xlab,Xname,Xt,Xfullprob] = load_processed_2DS_labels(dirname,t_str_start,t_str_stop,load_fullprob)

    if nargin < 4
        load_fullprob = 0;
    end

    fprintf('Load processed data...');

    if nargin > 1 && ~isempty(t_str_start) && ~isempty(t_str_stop)      
        
    else
        
        tmin = datenum('20000101000000','yyyymmddHHMMSS');
        tmax = datenum('20300101000000','yyyymmddHHMMSS');
        
    end

    file_list = dir(fullfile(dirname,'*.mat'));
    file_only_list = {file_list.name};
    file_list = fullfile(dirname,file_only_list);
    n_tot = numel(file_list);

    % loop over the snowflakes and save each feature in a vector of the main
    % structure data
    n_flakes = 0;
    X = zeros(n_tot,3);
    Xlab = cell(n_tot,1);
    Xname = cell(n_tot,1);
    Xt = zeros(n_tot,1);
    
    if load_fullprob
        load(file_list{1});
        tmp = get_struct_field(roi,'label_probs');
        n_class = numel(tmp);
        Xfullprob = zeros(n_tot,n_class);
    else
        Xfullprob = [];
    end
    
    i = 1;
    for k=1:length(file_list)

        load(file_list{k});
        % load only pictures within the time interval and large enough to process
        if (~isfield(roi,'tnum') || (roi.tnum >= tmin && roi.tnum <= tmax)) %&& ~roi.is2small

            j = 1;
            n_flakes = n_flakes + 1;
            
            Xname{i} = file_only_list{i};
            %disp(Xname{i});
            if isfield(roi,'tnum')
                Xt(i) = roi.tnum;
            else 
                Xt(i) = NaN;
            end

            X(i,j) = get_struct_field(roi,'label_ID'); Xlab{j} = 'label_ID'; j=j+1; %1
            X(i,j) = get_struct_field(roi,'area'); Xlab{j} = 'area'; j=j+1; %2
            X(i,j) = get_struct_field(roi,'Dmax'); Xlab{j} = 'Dmax'; j=j+1; %3
            X(i,j) = get_struct_field(roi,'label_ID2'); Xlab{j} = 'label_ID2'; j=j+1; %4
            X(i,j) = get_struct_field(roi,'roundness'); Xlab{j} = 'roundness'; j=j+1; %5
            X(i,j) = get_struct_field(roi,'compactness'); Xlab{j} = 'solidity (hull)'; j=j+1; %6

            if load_fullprob
                Xfullprob(i,:) = get_struct_field(roi,'label_probs');
            end
            
            i = i+1;

        end

    end
    
    
    
    fprintf('   Done! %u images found and processed. %u descriptors computed. \n',n_flakes,size(X,2));
    Xlab = Xlab';
    Xname = Xname';
    Xt = Xt';
    
end