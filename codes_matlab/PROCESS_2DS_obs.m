% === Core code to process airborne 2DS/HVPS probe images =================
%
% Process 2DS/HVPS probe images and create a matlab structure similar to
% what is done for MASC data. The code is designed to work with B&W binary
% silhouettes.
%
% Outputs are the cropped image and a structure containing process
% parameters and snowflake features.
%
% Author : Christophe Praz (christophe.praz@epfl.ch)
% Last update : May 2017
% =========================================================================
clear all;

tic;

% Load IO paths
%label.campaigndir = '/home/kiko/Documents/PhD/airlab/SampleImage/CPI/Aggregate';
label.campaigndir = '/home/praz/Documents/airlab/SampleImage/CPI/merged2';
label.outputdir_mat = '/home/praz/Documents/airlab/training_set/CPI2';
label.outputdir_img = '/home/praz/Documents/airlab/training_set/CPI2';
% ...

% Image processing parameters
process.probe = 'CPI';

process.min_hole_area = 10;
process.min_area_for_convex_hull = 15;
process.min_area_for_truncated_ellipse_fit = 40;
process.use_icpca = 1;
process.icpca_default = 0;
process.save_mat = 1;
process.save_img = 1;
process.illustration = 0;

% for 2DS (mostly)
process.smoothen_perim = 1;
% for CPI (mostly)
process.rm_white_boarders = true;
process.white_thresh = 0.95;
process.graythresh_multiplier_sup = 1.3;
process.graythresh_multiplier_inf1 = 0.3;
process.graythresh_multiplier_inf2 = 0.5;
process.graythresh_multiplier_inf3 = 0.7;
process.discard_noisy_img = true;
process.noise_threshold = 4.1; % 4.0-4.1 for adjusted img, 4.2 for raw img


% Upload particle images present in campaigndir
img_list = dir(fullfile(label.campaigndir,'*.png'));
img_list = {img_list.name};
img_list = img_list';

img_list2 = dir(fullfile(label.campaigndir,'*.jpg'));
img_list2 = {img_list2.name};
img_list2 = img_list2';

img_list = [img_list; img_list2];


fitted_OA_vec = [];
eq_radius = [];
n_skipped = 0;

for i=1:length(img_list)
    
    fprintf('Processing %s...\n',img_list{i});

    % filename
    DS.name = img_list{i}(1:end-4);
    DS.ext = img_list{i}(end-3:end);
    DS.probe = process.probe;
    
    % load data
    if strcmp(process.probe,'2DS') || strcmp(process.probe,'HVPS')
        tmp_mat = imread(fullfile(label.campaigndir,img_list{i})); % grayscale or truecolor
        DS.data = false(size(tmp_mat,1),size(tmp_mat,2));
        DS.data(tmp_mat(:,:,1)==0) = true;   
    
    % CPI probe requires a little bit more of processing as images are true colors and generally contain noise
    elseif strcmp(process.probe,'CPI')
        info = imfinfo(fullfile(label.campaigndir,img_list{i}));

        switch info.ColorType
            case 'grayscale'
                img_ini = im2double(imread(fullfile(label.campaigndir,img_list{i})));
                tmp_mat = img_ini;
            case 'truecolor'
                img_ini = imread(fullfile(label.campaigndir,img_list{i}));
                tmp_mat = rgb2gray(im2double(img_ini));
                [img_ini, map] = rgb2ind(img_ini,256);
            case 'indexed'
                [img_ini, map] = imread(fullfile(label.campaigndir,img_list{i}));
                tmp_mat = im2double(ind2gray(img_ini, map));
            otherwise
                error('invalid image type')
        end
        %tmp_mat = rgb2gray(im2double(imread(fullfile(label.campaigndir,img_list{i}))));
        tmp_mat_ini = tmp_mat;
        [h_ini,w_ini] = size(tmp_mat_ini);
        
    
        % detect and remove the "white boarders around the image"
        if process.rm_white_boarders
            med_col = median(tmp_mat,1);
            med_line = median(tmp_mat,2);
            idx_valid_col = find(med_col < process.white_thresh);
            idx_valid_line = find(med_line < process.white_thresh);
            tmp_mat = tmp_mat(idx_valid_line,idx_valid_col);
            
        end
    
        % alternative to remove background noise
        tmp_mat_corr = imadjust(tmp_mat);
        
        % for CPI only : detect "pure noise" images and discard them
        if process.discard_noisy_img
            
            img_HISE = fmeasure(tmp_mat_corr,'HISE',[]);
            
            if img_HISE < process.noise_threshold
                
                DS.status = 'noisy';
                
                
            end
    
        end
            
        
        % tricky part : binarize the image
        bw_mask_ini = [];
        bw_mask_ini{1} = ~imbinarize(tmp_mat_corr); % alternative 1 
        bw_mask_ini{2} = ~imbinarize(tmp_mat_corr,process.graythresh_multiplier_sup*graythresh(tmp_mat_corr)); % alternative 2
        bw_mask_ini{3} = ~imbinarize(tmp_mat_corr,process.graythresh_multiplier_inf1*graythresh(tmp_mat_corr)); % alternative 3
        bw_mask_ini{4} = ~imbinarize(tmp_mat_corr,process.graythresh_multiplier_inf2*graythresh(tmp_mat_corr)); % alternative 4
        bw_mask_ini{5} = ~imbinarize(tmp_mat_corr,process.graythresh_multiplier_inf3*graythresh(tmp_mat_corr)); % alternative 5
        %select best alternative based on particle shape complexity
        for j=1:numel(bw_mask_ini)
            tmp_roi = regionprops(bw_mask_ini{j},'Image','Area','Perimeter','PixelList');
            tmp_areas = [tmp_roi.Area];
            [area_max,idx_max] = max(tmp_areas);
            perim_max = tmp_roi(idx_max).Perimeter;
            score_alt(j) = perim_max/(2*sqrt(pi*area_max));
            perim = tmp_roi(idx_max).PixelList;

            idx_touch_edge = find(perim(:,1) == 1 | perim(:,1) == w_ini | perim(:,2) == 1 | perim == h_ini); 
            score_alt2(j) = 1 + numel(idx_touch_edge)/perim_max;
            score_alt3(j) = score_alt2(j) * score_alt(j); 
          
        end
            
        %score_alt_1 = sum(sum(bwperim(bw_mask_ini{1})))/(2*sqrt(pi*bwarea(bw_mask_ini{1})))
        %score_alt_2 = sum(sum(bwperim(bw_mask_ini{2})))/(2*sqrt(pi*bwarea(bw_mask_ini{2})))
        %score_alt_3 = sum(sum(bwperim(bw_mask_ini{3})))/(2*sqrt(pi*bwarea(bw_mask_ini{3})))
        %[~,idx_alt] = min([score_alt_1 score_alt_2 score_alt_3]);
        [~,idx_alt] = min(score_alt3);
        bw_mask_ini = bw_mask_ini{idx_alt};
        %bw_mask_ini = ~bw_mask_ini;
        %figure; imshow(bw_mask_ini);
        %bw_mask_fat = edge_detection(bw_mask_ini,2,0);
        bw_mask_ini_filled = imfill(bw_mask_ini,'holes');
        bw_holes_mask = logical(bw_mask_ini_filled - bw_mask_ini);
        bw_mask = bw_mask_ini;
        bw_mask(bw_holes_mask) = 0;
        %figure; imshow(bw_mask); 
        % at that point bw_mask should look like standard 2DS image
        DS.data = bw_mask;
        
    end

    % roi detection (normally there is just one main area, maybe 1-2 leftover pixels)
    all_roi = regionprops(DS.data,'Image','BoundingBox','SubarrayIdx', ...
         'PixelList','PixelIdxList','Perimeter','Area','MajorAxisLength', ...
         'Orientation','MinorAxisLength','Centroid'); 
     tmp_areas = [all_roi.Area];
     [max1,idx_selected] = max(tmp_areas);
     if max1 < process.min_area_for_convex_hull
         fprintf('Largest ROI found is smaller than %u pixels, particle classified as Small Particle (SP). \n',process.min_area_for_convex_hull);
         n_skipped = n_skipped + 1;
         DS.is2small = true;
     else
         DS.is2small = false;
     end
     % check the 2nd maximum to detect particles which might be "split"
     tmp_areas(idx_selected) = [];
     max2 = max(tmp_areas);

    mask = false(size(DS.data));
    mask(all_roi(idx_selected).PixelIdxList) = true;
    DS.data(~mask) = 0;
    B = bwboundaries(DS.data,'noholes');
    if numel(B) > 1
        disp('Warning: more than 1 particle detected on the image after the ROI selection procedure !!!');
    end
    B = B{1};
    DS.x_perim = B(:,2);
    DS.y_perim = B(:,1); 
    DS.perim = length(DS.x_perim);        

    % before cropping DS.data, determine if the flake is touching the boarder and if yes, by how much
    [h,w] = size(DS.data);
    idx = find(DS.x_perim == 1 | DS.x_perim == w | DS.y_perim == 1 | DS.y_perim == h); 
    idx2 = find(DS.x_perim == 2 | DS.x_perim == w-1 | DS.y_perim == 2 | DS.y_perim == h-1); 
    if ~isempty(idx)
        DS.touch_edge = 1;
        DS.touch_edge_Npix = length(idx);
        DS.touch_edge_ratio = length(idx)/DS.perim;
        DS.touch_edge_ratio_2 = length(idx2)/DS.perim;
    else
        DS.touch_edge = 0;
        DS.touch_edge_ratio = 0;
        DS.touch_edge_ratio_2 = 0;
    end

    % crop around the particle mask and clean the leftover pixels around the roi in the cropped rectangle
    crop_x = ceil(all_roi(idx_selected).BoundingBox(2));
    crop_y = ceil(all_roi(idx_selected).BoundingBox(1));
    crop_width = floor(all_roi(idx_selected).BoundingBox(3));
    crop_height = floor(all_roi(idx_selected).BoundingBox(4));
    DS.data = DS.data(crop_x:(crop_x+crop_height-1),crop_y:(crop_y+crop_width-1));  
    tmp_mat_cropped = tmp_mat(crop_x:(crop_x+crop_height-1),crop_y:(crop_y+crop_width-1)); 
    img_ini_cropped = img_ini(crop_x:(crop_x+crop_height-1),crop_y:(crop_y+crop_width-1)); 
    
    % compute the holes before smoothing (which requires filling)
    tmp_filled = imfill(DS.data,'holes');
    DS.bw_holes_mask = logical(tmp_filled - DS.data);
    
    
    if process.smoothen_perim
%         %figure;
        %subplot(121);
        %imshow(DS.bw_mask);
        % petit test lissage du perimetre
        DS.data = edge_detection(DS.data,3,0);
        tmp_rois = regionprops(DS.data,'Area','PixelIdxList');
        
        if numel(tmp_rois) > 1
            [~,idx_main] = max([tmp_rois.Area]);
            DS.data(:) = false;
            DS.data(tmp_rois(idx_main).PixelIdxList) = true;
        end
        % remove holes but keep perimeter
        tmp_perim = bwperim(DS.data);
        DS.data (DS.bw_holes_mask) = 0;
        DS.data(tmp_perim) = 1;
        
        %DS.bw_mask_filled = imfill(DS.bw_mask,'holes');
        %DS.data(DS.bw_holes_mask) = 0;

        %subplot(122);
        %imshow(DS.bw_mask);

    end

    % retrieve the BW mask of the selected particle (clean leftover pixels)
    % for binary images, the BW mask == data
    DS.bw_mask = DS.data;

    % retrieve horizontal/vertical dimensions from the selected particle
    DS.width = all_roi(idx_selected).BoundingBox(3);
    DS.height = all_roi(idx_selected).BoundingBox(4);
    DS.Dmean = 0.5*(DS.width + DS.height);

    % retrieve the fitted ellipse
    DS.E.a = all_roi(idx_selected).MajorAxisLength/2;
    DS.E.b = all_roi(idx_selected).MinorAxisLength/2;
    DS.E.theta = all_roi(idx_selected).Orientation;
    DS.centroid = all_roi(idx_selected).Centroid;

    % Generate mask of the snowflake without holes + compute N holes
    DS.bw_mask_filled = imfill(DS.bw_mask,'holes'); % mask filled     
    DS.bw_holes_mask = logical(DS.bw_mask_filled - DS.bw_mask);
    holes_roi = regionprops(DS.bw_holes_mask,'PixelIdxList','Perimeter','Area');
    holes_area = [holes_roi.Area];
    idx = find(holes_area >= process.min_hole_area);
    DS.nb_holes = numel(idx);
    
    if strcmp(process.probe,'CPI') 
        
        % Retrieve grayscale image of the cropped snowflake
        DS.data_gs = tmp_mat_cropped;
        DS.data_gs(~DS.bw_mask_filled) = 1; % 1 = white in grayscale
        % remapping and inversing colors in order to get something similar to MASC images
        DS.data_gs = uint8((1-DS.data_gs).*255);
        
        % Retrieve original "bluescale" image
        DS.data_bs = img_ini_cropped;
        DS.data_bs = ind2rgb(DS.data_bs,map);
        RGB_vals = [1 1 1];
        frame1 = DS.data_bs(:,:,1); frame1(~DS.bw_mask_filled) = RGB_vals(1);
        frame2 = DS.data_bs(:,:,2); frame2(~DS.bw_mask_filled) = RGB_vals(2);
        frame3 = DS.data_bs(:,:,3); frame3(~DS.bw_mask_filled) = RGB_vals(3);
        DS.data_bs(:,:,1) = frame1;
        DS.data_bs(:,:,2) = frame2;
        DS.data_bs(:,:,3) = frame3;

        % compute textural descriptors
        DS.mean_intens = mean(DS.data_gs(DS.bw_mask_filled))/255;
        DS.max_intens = double(max(DS.data_gs(DS.bw_mask_filled)))/255;
        range_array = rangefilt(DS.data_gs);
        DS.range_intens = mean(range_array(DS.bw_mask_filled))/255;
        DS.focus = DS.mean_intens * DS.range_intens;
        DS.std = std2(DS.data_gs(DS.bw_mask_filled));
        local_std = stdfilt(DS.data_gs);
        DS.local_std = mean(local_std(DS.bw_mask_filled)); 

        DS.lap = fmeasure(DS.data_gs,'LAPM',[]); 
        DS.hist_entropy = fmeasure(DS.data_gs,'HISE',[]);
        DS.wavs = fmeasure(DS.data_gs,'WAVS',[]);

        DS.H = haralick_props(DS.data_gs,0);
    
    end
    

    % the area + equivalent radius
    [DS.y DS.x] = find(DS.bw_mask_filled > 0);
    DS.area = sum(DS.bw_mask_filled(:));
    DS.area_porous = sum(DS.bw_mask(:));
    DS.n_true = sum(DS.data(:)); % number of 1 pixels
    DS.ntrue_ratio = DS.area_porous/DS.n_true;
    DS.eq_radius = sqrt(DS.area/pi);

    % recompute the perimeter on the cropped image (necessary to update the
    % pixel indices with respect to the new image dimensions)
    B = bwboundaries(DS.bw_mask_filled,'noholes');
    if numel(B) > 1
        disp('Warning: more than 1 particle detected on bw_mask_filled !!!');
    end
    B = B{1};
    DS.x_perim = B(:,2);
    DS.y_perim = B(:,1); 
    DS.perim = length(DS.x_perim);
    
    % the local centroid
    tmp = regionprops(DS.bw_mask_filled,'Centroid');
    DS.centroid_local = tmp.Centroid;
    DS.E.X0 = DS.centroid_local(1);
    DS.E.Y0 = DS.centroid_local(2); 
   
    if ~DS.is2small
        
        % the convex hull
        DS.hull = compute_convex_hull(DS.x,DS.y,DS.perim,0);
        
        if isnan(get_struct_field(DS.hull,'xh'))
            DS.C_out = NaN;
            DS.Dmax = NaN; % = max([DS.width, DS.height]);
            DS.D90 = NaN;
            DS.E_out = NaN;
            DS.Rect = NaN;
            DS.compactness = NaN;
            DS.roundness = NaN;
        else
            DS.C_out = fit_circle_around(DS.x_perim,DS.y_perim,DS.hull.xh,DS.hull.yh,0);
            [DS.Dmax,DS.Dmax_angle,DS.DmaxA,DS.DmaxB] = compute_Dmax(DS.hull.xh,DS.hull.yh,0); 
            DS.D90 = compute_D90(DS.bw_mask_filled,DS.Dmax,DS.Dmax_angle,0,0);
            DS.E_out = fit_ellipse_around(DS.x,DS.y,DS.hull.xh,DS.hull.yh,DS.E.theta*pi/180,0);
            DS.Rect = compute_rectangularity(DS.x,DS.y,DS.perim,0);
            DS.compactness = DS.area/(pi*DS.E_out.a*DS.E_out.b);
            DS.roundness = DS.area/DS.C_out.A;
        end
             
        % the complexity
        DS.complex = DS.perim/(2*pi*DS.eq_radius);
  
        % the inscribed ellipses
        DS.E_in = fit_ellipse_inside(DS.x,DS.y,DS.x_perim,DS.y_perim,DS.E.theta*pi/180,0);
        
        % the morphological skeleton
        DS.skel = skeleton_props(DS.data,0);
        
        % the fractal dimension
        DS.F = fractal_dim(DS.data,inf,0);
        DS.F_jac = 2*log(DS.perim/4)/log(DS.area); % [1,2]
        
        % the symmetry features
        DS.Sym = compute_symmetry_features(DS.bw_mask_filled,DS.Dmax,DS.eq_radius,0);
        
    else
        
        DS.hull = NaN;
        DS.complex = 1;
        DS.C_out = NaN;
        DS.Dmax = NaN; %max([DS.width, DS.height]);
        DS.Dmax_angle = NaN;
        DS.D90 = NaN;
        DS.E_in = NaN;
        DS.E_out = NaN;
        DS.Rect = NaN;
        DS.compactness = 1;
        DS.roundness = 1;
        DS.skel = NaN;
        DS.F = 2;
        DS.F_jac = 2;
        DS.Sym = NaN;
    end
  
    % the ellipse fitted to the particle outline without considering pixels touching the edges of the frame 
    % ---> this is used to compute a descriptor very helpful to detect truncated elliptical shapes, typically "Sphere" habit
    
    if DS.area >= process.min_area_for_truncated_ellipse_fit
    
        [h,w] = size(DS.bw_mask_filled);
        idx = find(DS.x_perim == 1 | DS.x_perim == w | DS.y_perim == 1 | DS.y_perim == h);
        tmp_xp = DS.x_perim; tmp_xp(idx) = [];
        tmp_yp = DS.y_perim; tmp_yp(idx) = [];
        if numel(tmp_xp) > 20
            try 
                DS.E2 = fit_ellipse_ls(tmp_xp,tmp_yp,w,h,0);
            catch
                DS.E2.BW_fit_area = NaN;
            end
        else
            DS.E2.BW_fit_area = NaN;
        end
    
    else
        
        DS.E2.BW_fit_area = NaN;
        
    end
    
    if ~isnan(DS.E2.BW_fit_area)
        dummy_BW = false(size(DS.bw_mask_filled)); %bw_mask_filled
        dummy_BW(DS.bw_mask_filled & DS.E2.BW_fit_area) = true;
        dummy_BW(~DS.bw_mask_filled & ~DS.E2.BW_fit_area) = true;
        DS.E2.fit_OA = sum(sum(dummy_BW))/(size(dummy_BW,1)*size(dummy_BW,2));
    else
        DS.E2.fit_OA = -1;
        fprintf('fit OA is : %2.2f \n',DS.E2.fit_OA);
        fprintf('Dmean / area are : %2.2f / %2.2f \n',DS.Dmean,DS.area);

    end
    
%     % additional descriptors computed with IC-PCA (Lindqvist et al. 2012)
%     masterPath = '/home/praz/Documents/IC-PCA/ref_setup/opetusmatsku/test_samples';
%     
% filePath = '/home/praz/Documents/IC-PCA/ref_setup/opetusmatsku/test_samples/Ros_4835_204560_2.png';
% inputDir = '/home/praz/Documents/IC-PCA/ref_setup/opetusmatsku/test_samples';
% outputDir = '/home/praz/Documents/IC-PCA/ref_setup/opetusmatsku/test_samples';
% filename = 'Ros_4835_204560_2.png';
% plots = true;
% debug = false; 
% enableFilters = false;
% imageExt = 'png';
% 
% 
% fileList = dir(fullfile(inputDir,['*' imageExt])); 
% paramfile = fullfile(outputDir,'class_parameters.out'); % file containing the classification parameters
% namesfile = fullfile(outputDir,'filenames.out');
% aste = 1;

    icpca = [];

    if process.use_icpca
        
        if DS.is2small
            
            DS.icpca.a = NaN;
            DS.icpca.b = NaN;
            
        else
        
        icpca.filePath = fullfile(label.campaigndir,strcat(DS.name,DS.ext));
        icpca.outputDir = '/home/praz/Documents/IC-PCA/figs/CPI';
        icpca.plots = true;
        icpca.debug = false;
        icpca.enableFilters = false;
        if process.icpca_default % 1:use raw image as input | 0:use DS.data as input 
            [icpca.x,icpca.y,icpca.amax,icpca.area_no_holes,icpca.name,icpca.cout,icpca.centroid,icpca.adiff,icpca.imsize,icpca.edgepx,icpca.xtremepx] = get_perimeter(icpca.filePath, icpca.outputDir, icpca.plots, icpca.debug, icpca.enableFilters);
        else
            [icpca.x,icpca.y,icpca.amax,icpca.area_no_holes,icpca.name,icpca.cout,icpca.centroid,icpca.adiff,icpca.imsize,icpca.edgepx,icpca.xtremepx] = get_perimeter(icpca.filePath, icpca.outputDir, icpca.plots, icpca.debug, icpca.enableFilters,double(DS.data));
        end
        DS.icpca = get_parameters(DS.name, label.campaigndir, icpca.outputDir, icpca.outputDir, icpca.plots, 1, icpca.x, icpca.y, icpca.amax, icpca.area_no_holes, icpca.name, icpca.cout, icpca.centroid, icpca.adiff, icpca.imsize, icpca.edgepx, icpca.xtremepx, icpca.enableFilters, icpca.outputDir);
        
        end
        
    end
    
    if process.illustration
        
        figure;
        subplot(221); imshow(tmp_mat_ini); xlabel('Raw image');
        subplot(222); imshow(DS.data_gs); xlabel('Final image');
        subplot(223); imshow(DS.data); xlabel('BW mask');

    end
    
    % save results in a .mat structure and the cropped particle in a .png
    if process.save_mat

        % the .mat structure
        if ~exist(label.outputdir_mat,'dir')
            mkdir(label.outputdir_mat);
        end
        save_mat(label.outputdir_mat,strcat(DS.name,'.mat'),DS);
        
    end
    
    if process.save_img

        % the .png cropped particle image
        if ~exist(label.outputdir_img,'dir')
            mkdir(label.outputdir_img);
        end
        if strcmp(process.probe,'2DS') || strcmp(process.probe,'HVPS')
            img2save = DS.data;
            img2save(DS.data) = false; % reverse colors
            img2save(~DS.data) = true;
            imwrite(img2save.*255,fullfile(label.outputdir_img,strcat(DS.name,'.png')),'png','BitDepth', 8);
            
        elseif strcmp(process.probe,'CPI') 
            img2save = DS.data_bs;
            imwrite(img2save,fullfile(label.outputdir_img,strcat(DS.name,'.png')),'png','BitDepth',8);
            
        end

    end

    % delete the saved structure
    if i < numel(img_list) 
        clear DS;
    end
   
end

t1 = toc;
fprintf('Done ! It took %2.2f secs to process %u images. \n',t1,length(img_list));
    

%% illustrations for debugging 
if 0
    figure; hold on; grid on; box on;
    histogram(fitted_OA_vec(1:1008),[0.5:0.005:1]);
    histogram(fitted_OA_vec(1009:end),[0.5:0.005:1]);
    xlabel('OA ratio truncated ellipse/particle area');
    ylabel('count');
    legend('Others','Spheres');
    set(gca,'Fontsize',14);

    figure; hold on; grid on; box on;
    plot(eq_radius(1:1008),fitted_OA_vec(1:1008),'ko','MarkerFaceColor','b','MarkerSize',3);
    plot(eq_radius(1009:end),fitted_OA_vec(1009:end),'ko','MarkerFaceColor','r','MarkerSize',3);
end


