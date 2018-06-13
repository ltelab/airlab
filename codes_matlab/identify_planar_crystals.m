% small piece of code to isolate potential planar crystals candidates
clearvars;
data_dir = '/home/praz/Documents/2DS_data/2DSinputDir_png_20151112_194000-195000/Pla_proc';
target_dir = '/home/praz/Documents/2DS_data/2DSinputDir_png_20151112_194000-195000/Pla_proc_filtered';
load('potential_planar_crystals_2.mat');
n_max = 54;

for i=1:n_max
    
    fprintf('%u/%u\n',i,n_max);
    current_name = Xname_hope_cool{i};
    current_name = strcat(current_name(1:end-4),'.png');
    img_name = dir(fullfile(data_dir,'**',current_name));
    
    current_name_mat = strcat(current_name(1:end-4),'.mat');
    
    if numel(img_name) ~= 1
        fprintf('Warning : particle not found, skipped \n');
        continue;
        
    else
        img_fullpath = fullfile(img_name.folder,img_name.name);
        copyfile(img_fullpath,target_dir);
        
    end
        
   
    
end