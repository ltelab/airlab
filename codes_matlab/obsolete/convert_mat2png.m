% small script to save .mat images to .png

clear all; close all;

data_dir = '../Saisai_2DS_classification_case/inputDir';

img_list = dir(fullfile(data_dir,'*.mat'));
img_list = {img_list.name};
img_list = img_list';

for i=1:length(img_list)
    load(fullfile(data_dir,img_list{i}));
    name = img_list{i};
    name = name(1:end-4);
    name = strcat(name,'.png');
    imwrite(particle_mask.*255,fullfile(data_dir,name),'png','BitDepth',8);
end
