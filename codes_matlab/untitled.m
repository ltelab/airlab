% small piece of code to illustrate the concept of invariant angle and line
% segment (IC-PCA) on different shapes

clearvars; close all;

            tmp_mat = imread(fullfile(img_path,img_current)); % grayscale or truecolor
            DS.data = false(size(tmp_mat,1),size(tmp_mat,2));
            DS.data(tmp_mat(:,:,1)==0) = true; 
            DS.raw_data = DS.data;
            [DS.h_ini,DS.w_ini] = size(DS.raw_data);