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
clearvars;

tic;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%% USER PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% path to the raw images
label.campaigndir = '../SampleImage/CPI/merged_all';

% path to the directory to save processed matfiles (.mat)
label.outputdir_mat = '../training_set/CPI_smooth0_icpca1';

% path to the directory to save processed images (.png)
label.outputdir_img = '../training_set/CPI_smooth0_icpca1';

% imaging probe. Can be : 2DS / HVPS / CPI
process.probe = 'CPI';

% input format : 'image' (.png) or 'matfile' (.mat) (matfile supported only for 2DS probe
process.input_img_type = 'image';

% run the code in parallel (true) or not (false)
process.parallel = true;

% save processed matfiles ? (true/false)
process.save_mat = true;

% save processed images ? (true/false)
process.save_img = true;

% illustration of the image cropping procedure ? (true/false, used mostly for debugging)
process.illustration = 0;

% some other stuff
process.min_hole_area = 10; % 10 default
process.min_area_for_convex_hull = 15; % 15 default
process.min_area_for_truncated_ellipse_fit = NaN; % feature currently disabled
process.use_icpca = true; % use ICPCA code from Lindquvist et al.

% some technical parameters used to identify background threshold in CPI images
process.rm_white_boarders = false;
process.white_thresh = 0.95;
process.graythresh_multiplier_sup = 1.3;
process.graythresh_multiplier_inf1 = 0.7;
process.graythresh_multiplier_inf2 = 0.5;
process.graythresh_multiplier_inf3 = 0.3;
process.noise_threshold = 4.1; % 4.0-4.1 for adjusted img, 4.2 for raw img

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if strcmp(process.probe,'2DS') || strcmp(process.probe,'HVPS')
    process.icpca_default = true;
    process.smoothen_perim = false;
    process.discard_noisy_img = false;

elseif strcmp(process.probe,'CPI')
    process.icpca_default = 1; %default = false
    process.smoothen_perim = 0; % default = true
    process.discard_noisy_img = false;
    
end

% Upload particle images present in campaigndir
if strcmp(process.input_img_type,'image')
    img_list = dir(fullfile(label.campaigndir,'*.png'));
    img_list = {img_list.name};
    img_list = img_list';
    img_list2 = dir(fullfile(label.campaigndir,'*.jpg'));
    img_list2 = {img_list2.name};
    img_list2 = img_list2';
    img_list = [img_list; img_list2];
    
elseif strcmp(process.input_img_type,'matfile')
    img_list = dir(fullfile(label.campaigndir,'*.mat'));
    img_list = {img_list.name};
    img_list = img_list';
    
end
 
n_tot = numel(img_list);

if process.parallel
    % if you have less than 16 threads matlab will take as many as possible
    n_workers = 16; 
    fprintf('The code will run on parallel on your computer, the progress bar cannot be displayed\n');
else
    n_workers = 0;
end

parfor (i=1:n_tot,n_workers) 
    
    if n_workers == 0
        fprintf('%u/%u : %s\n',i,n_tot,img_list{i});
    end
    process_OAP_image(img_list{i},label,process);
    
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


