% first attempt to detect pure noise images in CPI data
clearvars; close all;

noisy_dir = '/home/praz/Documents/airlab/SampleImage/CPI/noise';
good_dir = '/home/praz/Documents/airlab/SampleImage/CPI/merged1';

HISE_thresh = 4.1;

all_dirs = {noisy_dir; good_dir};

for k=1:numel(all_dirs)

    img_list = dir(fullfile(all_dirs{k},'*.png'));
    img_list = {img_list.name};
    img_list = img_list';
    img_list2 = dir(fullfile(all_dirs{k},'*.jpg'));
    img_list2 = {img_list2.name};
    img_list2 = img_list2';
    img_list = [img_list; img_list2];

    white_thresh = 0.95;

    for i=1:numel(img_list)

        info = imfinfo(fullfile(all_dirs{k},img_list{i}));

        switch info.ColorType
            case 'grayscale'
                img_ini = im2double(imread(fullfile(all_dirs{k},img_list{i})));
                img = img_ini;
            case 'truecolor'
                img_ini = imread(fullfile(all_dirs{k},img_list{i}));
                img = rgb2gray(im2double(img_ini));
                [img_ini, map] = rgb2ind(img_ini,256);
            case 'indexed'
                [img_ini, map] = imread(fullfile(all_dirs{k},img_list{i}));
                img = im2double(ind2gray(img_ini, map));
            otherwise
                error('invalid image type')
        end

        % detect and remove the "white boarders around the image" (probably not necessary)
        med_col = median(img,1);
        med_line = median(img,2);
        idx_valid_col = find(med_col < white_thresh);
        idx_valid_line = find(med_line < white_thresh);
        img = img(idx_valid_line,idx_valid_col);
        
        img = imadjust(img);
        
%         try
%             img = adapthisteq(img);
%         catch
%             disp('lol');
%             end

        % noise descriptor 1 (gaussian filtering)
        %tmp_mean = mean(size(img{i}));
        sigma_blur = 0.25;
        img_blur = imgaussfilt(img,sigma_blur);

        diff = img-img_blur;
        stdv{k}(i,1) = std2(diff);
        meanv{k}(i,1) = mean(diff(:));
        meanv_abs{k}(i,1) = mean(abs(diff(:)));

        
        local_std{k}(i,1) = mean(mean(stdfilt(img)));
        
        
        HISE{k}(i,1) = fmeasure(img,'HISE',[]);
        LAPM{k}(i,1) = fmeasure(img,'LAPM',[]);
        WAVS{k}(i,1) = fmeasure(img,'WAVS',[]);
        
        H{k}(i,1) = haralick_props(img);
       
        
        %figure;
        %subplot(131); imshow(img); title('Raw');
        %subplot(132); imshow(img_blur); title('Filtered');
        %subplot(133); imshow(diff); title('Diff');

    end
    
end

%%

figure; 
subplot(221); hold on; box on;
title('mean bias');
histogram(meanv{1});
histogram(meanv{2});
legend('noise','good');

subplot(222); hold on; box on;
title('diff std');
histogram(stdv{1});
histogram(stdv{2});
legend('noise','good');

subplot(223); hold on; box on;
title('mean abs bias');
histogram(meanv_abs{1});
histogram(meanv_abs{2});
legend('noise','good');

subplot(224); hold on; box on;
title('mean local std');
histogram(local_std{1});
histogram(local_std{2});
legend('noise','good');

figure; 
subplot(131); hold on; box on;
title('histogram entropy');
histogram(HISE{1});
histogram(HISE{2});
plot([HISE_thresh HISE_thresh],[0 350],'r--');
legend('noise','good');

subplot(132); hold on; box on;
title('Laplacian energy');
histogram(LAPM{1});
histogram(LAPM{2});
legend('noise','good');

subplot(133); hold on; box on;
title('WAVS');
histogram(WAVS{1});
histogram(WAVS{2});
legend('noise','good');

% Haralick
figure; 
subplot(221); hold on; box on;
title('H Contrast');
histogram([H{1}.Contrast],'Normalization','pdf');
histogram([H{2}.Contrast],'Normalization','pdf');
legend('noise','good');

subplot(222); hold on; box on;
title('H Energy');
histogram([H{1}.Energy]);
histogram([H{2}.Energy]);
legend('noise','good');

subplot(223); hold on; box on;
title('H Homogeneity');
histogram([H{1}.Homogeneity]);
histogram([H{2}.Homogeneity]);
legend('noise','good');

subplot(224); hold on; box on;
title('H H Correlation');
histogram([H{1}.Correlation]);
histogram([H{2}.Correlation]);
legend('noise','good');









