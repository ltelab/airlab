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
label.campaigndir = '/home/praz/Documents/CPI_data/tmp1';
label.outputdir_mat = '/home/praz/Documents/CPI_data/tmp2';
label.outputdir_img = '/home/praz/Documents/CPI_data/tmp2';
% ...

% Image processing parameters
process.probe = 'CPI';

process.min_hole_area = 1;
process.min_area_for_convex_hull = 15;
process.min_area_for_truncated_ellipse_fit = 40;
process.use_icpca = 1;
process.icpca_default = 0;
process.save_results = 1;
process.illustration = 0;

% for 2DS (mostly)
process.smoothen_perim = 0;
% for CPI (mostly)
process.rm_white_boarders = false;
process.white_thresh = 0.95;
process.graythresh_multiplier = 1.3;


% Upload particle images present in campaigndir
img_list = dir(fullfile(label.campaigndir,'*.png'));
img_list = {img_list.name};
img_list = img_list';

fitted_OA_vec = [];
eq_radius = [];
n_skipped = 0;

for i=1:length(img_list)
    
    fprintf('Processing %s...\n',img_list{i});

    % filename
    DS.name = img_list{i};
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
                tmp_mat = imread(fullfile(label.campaigndir,img_list{i}));
            case 'truecolor'
                tmp_mat = rgb2gray(im2double(imread(fullfile(label.campaigndir,img_list{i}))));
            case 'indexed'
                [X, map] = imread(fullfile(label.campaigndir,img_list{i}));
                tmp_mat = im2double(ind2gray(X, map));
            otherwise
                error('invalid image type')
        end
        %tmp_mat = rgb2gray(im2double(imread(fullfile(label.campaigndir,img_list{i}))));
        tmp_mat_ini = tmp_mat;
    
        % detect and remove the "white boarders around the image"
        if process.rm_white_boarders
            med_col = median(tmp_mat,1);
            med_line = median(tmp_mat,2);
            idx_valid_col = find(med_col < process.white_thresh);
            idx_valid_line = find(med_line < process.white_thresh);
            tmp_mat = tmp_mat(idx_valid_line,idx_valid_col);
            
        end
    
        %figure; imshow(tmp_mat);

        % optional smoothing (bad)
        %tmp_mat = conv2(tmp_mat,double(ones(3)/9),'full'); % by default : 'valid' 
        % figure; imshow(tmp_mat);
    
        % tricky part : binarize the image
        bw_mask_ini = [];
        bw_mask_ini{1} = imfill(~imbinarize(tmp_mat,graythresh(tmp_mat)),'holes'); % alternative 1 
        bw_mask_ini{2} = imfill(~imbinarize(tmp_mat,process.graythresh_multiplier*graythresh(tmp_mat)),'holes'); % alternative 2
        %select best alternative based on particle shape complexity
        score_alt_1 = sum(sum(bwperim(bw_mask_ini{1})))/(2*sqrt(pi*bwarea(bw_mask_ini{1})));
        score_alt_2 = sum(sum(bwperim(bw_mask_ini{2})))/(2*sqrt(pi*bwarea(bw_mask_ini{2})));
        [~,idx_alt] = min([score_alt_1 score_alt_2]);
        bw_mask_ini = bw_mask_ini{idx_alt};
        %bw_mask_ini = ~bw_mask_ini;
        %figure; imshow(bw_mask_ini);
        bw_mask_fat = edge_detection(bw_mask_ini,10,0);
        bw_mask_ini_filled = imfill(bw_mask_ini,'holes');
        bw_holes_mask = logical(bw_mask_ini_filled - bw_mask_ini);
        bw_mask = bw_mask_fat;
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

     % here is a small code to toroughly retrieve the particle outline
     % using an edge dilation/erosion procedure. This feature is currently
     % disabled as the classification is performing better without 
%      if max2 > 0.1*max1
%          fprintf('Warning: second largest ROI detected measures %2.1f of the main ROI area. \n',max2/max1*100);
% 
%         % test (optimal ?)
%         DS.data_original = DS.data;
%         DS.data_bulk = edge_detection(DS.data);
%         % detect the biggest roi in data_bulk
%         all_roi = regionprops(DS.data_bulk,'Image','BoundingBox','SubarrayIdx', ...
%                 'PixelList','PixelIdxList','Perimeter','Area','MajorAxisLength', ...
%                 'Orientation','MinorAxisLength','Centroid'); 
%         tmp_areas = [all_roi.Area];
%         [max1,idx_selected] = max(tmp_areas);
%         % keep only the biggest roi
%         DS.data_bulk = false(size(DS.data_bulk));
%         DS.data_bulk(all_roi(idx_selected).PixelIdxList) = true;
%         % retrieve the perimeter of this roi
%         B = bwboundaries(DS.data_bulk,'noholes');
%         if numel(B) > 1
%             disp('Warning: more than 1 particle detected on data_bulk !!!');
%         end
%         B = B{1};
%         DS.x_perim = B(:,2);
%         DS.y_perim = B(:,1); 
%         DS.perim = length(DS.x_perim);
%         DS.data_bulk_perim = bwperim(DS.data_bulk);
%         % recreate the flake using this perimeter
%         DS.data = DS.data_bulk;
%         DS.data(~DS.data_original) = false;
%         DS.data(DS.data_bulk_perim) = true;
% 
%     end

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
    if ~isempty(idx)
        DS.touch_edge = 1;
        DS.touch_edge_Npix = length(idx);
        DS.touch_edge_ratio = length(idx)/DS.perim;
    else
        DS.touch_edge = 0;
        DS.touch_edge_ratio = 0;
    end

    % crop around the particle mask and clean the leftover pixels around the roi in the cropped rectangle
    crop_x = ceil(all_roi(idx_selected).BoundingBox(2));
    crop_y = ceil(all_roi(idx_selected).BoundingBox(1));
    crop_width = floor(all_roi(idx_selected).BoundingBox(3));
    crop_height = floor(all_roi(idx_selected).BoundingBox(4));
    DS.data = DS.data(crop_x:(crop_x+crop_height-1),crop_y:(crop_y+crop_width-1));  
    tmp_mat_cropped = tmp_mat(crop_x:(crop_x+crop_height-1),crop_y:(crop_y+crop_width-1)); 

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
    
    % Update grayscale image of the snowflake
    DS.data_gs = tmp_mat_cropped;
    DS.data_gs(~DS.bw_mask_filled) = 1; % 1 = white in grayscale
    
    
    
    if process.smoothen_perim
        %figure;
        %subplot(121);
        %imshow(DS.bw_mask);

        % petit test lissage du perimetre
        DS.data = edge_detection(DS.data,0);
        tmp_rois = regionprops(DS.data,'Area','PixelIdxList');
        if numel(tmp_rois) > 1
            [~,idx_main] = max([tmp_rois.Area]);
            DS.data(:) = false;
            DS.data(tmp_rois(idx_main).PixelIdxList) = true;
        end
        DS.bw_mask_filled = imfill(DS.bw_mask,'holes');
        DS.data(DS.bw_holes_mask) = 0;

        %subplot(122);
        %imshow(DS.bw_mask);

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
        
        icpca.filePath = fullfile(label.campaigndir,DS.name);
        icpca.outputDir = '/home/praz/Documents/IC-PCA/figs/2DS';
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
    if process.save_results

        % the .mat structure
        if ~exist(label.outputdir_mat,'dir')
            mkdir(label.outputdir_mat);
        end
        save_mat(label.outputdir_mat,strcat(DS.name(1:end-4),'.mat'),DS);

        % the .png cropped particle image
        if ~exist(label.outputdir_img,'dir')
            mkdir(label.outputdir_img);
        end
        if strcmp(process.probe,'2DS') || strcmp(process.probe,'HVPS')
            img2save = DS.data;
            img2save(DS.data) = false; % reverse colors
            img2save(~DS.data) = true;
            imwrite(img2save.*255,fullfile(label.outputdir_img,DS.name),'png','BitDepth', 8);
            
        elseif strcmp(process.probe,'CPI')
            img2save = DS.data_gs;
            imwrite(img2save,fullfile(label.outputdir_img,DS.name),'png','BitDepth', 8);
            
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


