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

% Load IO paths
%label.campaigndir = '/home/kiko/Documents/PhD/airlab/SampleImage/CPI/Aggregate';
label.campaigndir = '/home/kiko/Documents/PhD/airlab/SampleImage/2DS/merged_all';
label.outputdir_mat = '/home/kiko/Documents/PhD/airlab/training_set/2DS_test';
label.outputdir_img = '/home/kiko/Documents/PhD/airlab/training_set/2DS_test';
% ...

% Image processing parameters
process.probe = '2DS';

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
process.rm_white_boarders = false;
process.white_thresh = 0.95;
process.graythresh_multiplier_sup = 1.3;
process.graythresh_multiplier_inf1 = 0.7;
process.graythresh_multiplier_inf2 = 0.5;
process.graythresh_multiplier_inf3 = 0.3;
process.discard_noisy_img = false;
process.noise_threshold = 4.1; % 4.0-4.1 for adjusted img, 4.2 for raw img


% Upload particle images present in campaigndir
img_list = dir(fullfile(label.campaigndir,'*.png'));
img_list = {img_list.name};
img_list = img_list';

img_list2 = dir(fullfile(label.campaigndir,'*.jpg'));
img_list2 = {img_list2.name};
img_list2 = img_list2';

img_list = [img_list; img_list2];
n_tot = numel(img_list);

parfor i=1:n_tot
    
    fprintf('%u : %s\n',n_tot,img_list{i});
    process_OAP_image(img_list{i},label,process); 
    %n_current = 1;
   
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


