% small script to assign labels to airborne probes particles
clear; close all;

dirname = '/ltedata/MASC/OAP/airlab/training_set/subclassification/bulrosagg_vs_others/all';

probe = 'CPI';
save_results = true;

file_list = dir(fullfile(dirname,'*.mat'));
file_only_list = {file_list.name};
file_list = fullfile(dirname,file_only_list);

for i=1:length(file_list)
    
    load(file_list{i});
    prefix = file_only_list{i}(1:3);
    %prefix = roi.name(1:3);
    
    switch(prefix)
        case 'Agg'
            roi.label_name = 'Aggregate';
            roi.label_ID = 1;
        case 'Col'
            roi.label_name = 'Column';
            roi.label_ID = 2;
        case 'Gra'
            roi.label_name = 'Graupel';
            roi.label_ID = 3;
        case 'Ros'
            roi.label_name = 'Rosette';
            roi.label_ID = 4;
        case 'Sph'
            roi.label_name = 'Sphere';
            roi.label_ID = 5;
        case 'Oth'
            roi.label_name = 'Other';
            roi.label_ID = 6;
        case 'Pla'
            roi.label_name = 'Plate';
            if strcmp(probe,'CPI')
                roi.label_ID = 6;
            elseif strcmp(probe,'2DS')
                roi.label_ID = 7;
            end
        % subclassification training set
        case 'BR_'
            roi.label_name = 'Aggregate';
            roi.label_ID = 1;            
            roi.label_name2 = 'Aggregate of bullet rosettes';
            roi.label_ID2 = 0;
        case 'O_A'
            roi.label_name = 'Aggregate';
            roi.label_ID = 1;
            roi.label_name2 = 'Other aggregate';
            roi.label_ID2 = 1;   
        otherwise
            roi.label_name = 'Unknown/Error';
            roi.label_ID = -1;
            fprintf('Warning : unexpected name for file %s \n',file_list{i});
    end
    
    if save_results
        save(file_list{i},'roi');
    end
    
end