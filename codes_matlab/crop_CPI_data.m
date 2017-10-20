% crop_CPI
clear all; close all;

datadir = '/home/praz/Documents/CPI_data/20151112/1854';
outdir = '/home/praz/Documents/CPI_data/tmp1';

all = dir(fullfile(datadir,'**','_CPI__20151112_192955_964_.png'));
all_name = {all.name}';
all_folder = {all.folder}';

for k = 1:numel(all_name)

    fprintf('Processing composite %s... ',all_name{k});
    filename = fullfile(all_folder{k},all_name{k});
    info = imfinfo(filename);

    switch info.ColorType
        case 'grayscale'
            im = imread(filename);
        case 'truecolor'
            im = rgb2gray(im2double(filename));
        case 'indexed'
            [X, map] = imread(filename);
            im = im2double(ind2gray(X, map));
        otherwise
            error('invalid image type')
    end



    BW_im = 1-im;

    BW_im = imfill(imbinarize(BW_im,0.01),'holes');

    se = strel('square',10);

    BW_im_dilated = imopen(BW_im,se);

    ROI = regionprops(BW_im_dilated,'BoundingBox');

    for i=1:numel(ROI)

        crop_x = ceil(ROI(i).BoundingBox(2));
        crop_y = ceil(ROI(i).BoundingBox(1));
        crop_width = floor(ROI(i).BoundingBox(3));
        crop_height = floor(ROI(i).BoundingBox(4));
        sub_roi = im(crop_x:(crop_x+crop_height-1),crop_y:(crop_y+crop_width-1));
        RGB = ind2rgb(sub_roi.*255,map);
        imwrite(sub_roi,fullfile(outdir,sprintf('part%u_%s',i,all_name{k})));

        %sub_roi = im(ROI(i).PixelIdxList);
        %figure;imshow(RGB,map);    

    end
    
    fprintf('Done !\n');
    
end


% %roi = regionprops(im);
% 
%     se0 = strel('line', strel_size, 0); %horz floor(1.5*process.linefill)/flake.imfinfo.XRes
%     se90 = strel('line', strel_size, 90); %vert % 9.79 floor(1.5*process.linefill)/flake.imfinfo.XRes
%     edged = edge(extended,'Sobel',0.008); %'log' is better but significantly slower       
% 
%     
%     % Dilate and erode in order to blurry internal complexity
%     dilated = imdilate(edged, [se0 se90]);
%     filled = imfill(dilated,'holes');


% figure; subplot(131); imshow(im); subplot(132); imshow(BW_im); subplot(133); imshow(BW_im_dilated);
