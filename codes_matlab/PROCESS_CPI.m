% === Core code to process airborne CPI probe images =================
%
% Process CPI probe images and create a matlab structure similar to
% what is done for MASC data. The code is designed to work with RGB .png
% CPI images
%
% Outputs are the cropped image and a structure containing process
% parameters and snowflake features.
%
% Author : Christophe Praz (christophe.praz@epfl.ch)
% Last update : October 2017
% =========================================================================
clear all; close all;

%tic;

% Load IO paths
label.campaigndir = '/home/praz/Documents/airlab/SampleImage/CPI/merged';
label.outputdir_mat = '/home/praz/Documents/airlab/training_set/CPI';
label.outputdir_img = '/home/praz/Documents/airlab/training_set/CPI';
% ...

% Image processing parameters
process.white_thresh = 0.95;
process.min_hole_area = 1;
process.min_area_for_convex_hull = 15;
process.min_area_for_truncated_ellipse_fit = 40;
process.save_results = 1;
process.probe = 'CPI';
process.use_icpca = 1;
process.icpca_default = 0;
process.smoothen_perim = 1;
% ...

% Upload particle images present in campaigndir
img_list = dir(fullfile(label.campaigndir,'Agg_28.png'));
img_list = {img_list.name};
img_list = img_list';

fitted_OA_vec = [];
eq_radius = [];
n_skipped = 0;

for i=1:length(img_list)
    
    fprintf('Processing %s...\n',img_list{i});

    roi.name = img_list{i};
    roi.probe = process.probe;
    
    tmp_mat = rgb2gray(im2double(imread(fullfile(label.campaigndir,img_list{i}))));
    
    figure; imshow(tmp_mat);
    
    % detect and remove the "white boarders around the image"
    med_col = median(tmp_mat,1);
    med_line = median(tmp_mat,2);
    idx_valid_col = find(med_col < process.white_thresh);
    idx_valid_line = find(med_line < process.white_thresh);
    tmp_mat = tmp_mat(idx_valid_line,idx_valid_col);
    
    figure; imshow(tmp_mat);

    % image smoothing (?)
    %tmp_mat = conv2(tmp_mat,double(ones(3)/9),'full'); % by default : 'valid' 
    % figure; imshow(tmp_mat);
    
    roi.bw_mask = imbinarize(tmp_mat,1.3*graythresh(tmp_mat));
    roi.bw_mask = ~roi.bw_mask;
    figure; imshow(roi.bw_mask);
    roi.bw_mask_filled = edge_detection(roi.bw_mask,1);
    tmp_filled = imfill(roi.bw_mask,'holes');
    roi.bw_holes_mask = logical(tmp_filled - roi.bw_mask);
    roi.bw_mask = roi.bw_mask_filled;
    roi.bw_mask(roi.bw_holes_mask) = 0;
    figure; imshow(roi.bw_mask);
    
    % roi detection (normally there is just one main area, maybe 1-2 leftover pixels)
    all_roi = regionprops(roi.bw_mask,'Image','BoundingBox','SubarrayIdx', ...
         'PixelList','PixelIdxList','Perimeter','Area','MajorAxisLength', ...
         'Orientation','MinorAxisLength','Centroid'); 
     tmp_areas = [all_roi.Area];
     [max1,idx_selected] = max(tmp_areas);
     if max1 < process.min_area_for_convex_hull
         fprintf('Largest ROI found is smaller than %u pixels, particle classified as Small Particle (SP). \n',process.min_area_for_convex_hull);
         n_skipped = n_skipped + 1;
         roi.is2small = true;
     else
         roi.is2small = false;
     end

    
    
    %    DS.bw_mask_filled = imfill(DS.bw_mask,'holes'); % mask filled     
    %DS.bw_holes_mask = logical(DS.bw_mask_filled - DS.bw_mask);
    %holes_roi = regionprops(DS.bw_holes_mask,'PixelIdxList','Perimeter','Area');
    %holes_area = [holes_roi.Area];
    %idx = find(holes_area >= process.min_hole_area);
    %DS.nb_holes = numel(idx);

% % Luodaan maski käyttäen thresholdia ja sopivaa kerrointa (inputparametri).
% mask=im2bw(im,graythresh(im)*threshcoef);
    
    
    % filename
%     DS.name = img_list{i};
%     DS.probe = process.probe;
%     % data (bw)
%     tmp_mat = imread(fullfile(label.campaigndir,img_list{i}));
%     DS.data = false(size(tmp_mat,1),size(tmp_mat,2));
%     DS.data(tmp_mat(:,:,1)==0) = true;
%     DS.n_true = sum(DS.data(:)); % number of 1 pixels 
% 



end