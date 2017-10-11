% small script to assign labels to airborne probes particles
clear; close all;

dirname = '/home/praz/Documents/airlab/training_set/2DS_smooth';

save_results = true;

file_list = dir(fullfile(dirname,'*.mat'));
file_only_list = {file_list.name};
file_list = fullfile(dirname,file_only_list);

for i=1:length(file_list)
    
    load(file_list{i});
    
    prefix = roi.name(1:3);
    
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
        otherwise
            roi.label_name = 'Unknown/Error';
            roi.label_ID = -1;
            fprintf('Warning : unexpected name for file %s \n',file_list{i});
    end
    
    if save_results
        save(file_list{i},'roi');
    end
    
end