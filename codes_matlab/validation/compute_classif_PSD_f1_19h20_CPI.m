clearvars; close all;

datadir = '/ltedata/MASC/OAP/airlab/training_set/CPI_smooth0_icpca0'; % CPI
savedir = '/ltedata/MASC/OAP/airlab/figures/bright_spot/training_gra';
[X, Xlab, Xname, Xt, Xfullprob] = load_processed_2DS_labels(datadir,[],[],true);
y = X(:,1);%load_2DS_labels(datadir{i});
Nini = numel(y);
Dmax_tmp = X(:,3);
area_tmp = X(:,2);
round = X(:,5);
compact = X(:,6); compact(compact > 1) = 1;
probMAX = max(Xfullprob,[],2);
% Dmax estimated from area for particles too small
Dmax_tmp(isnan(Dmax_tmp) & ~isnan(area_tmp)) = sqrt(area_tmp(isnan(Dmax_tmp) & ~isnan(area_tmp)));
fprintf('\n*** %s cloud probe data loaded *****\n','CPI');
fprintf('Number of images found : %u \n',Nini);
fprintf('************************************\n\n');    

%%
k=0;
probTHRESH = 0.;
idx = find(round>0.8 & Dmax_tmp>10 & y > 0);%find(y == 5 & probMAX > probTHRESH);
Xname = Xname(idx);
for i=1:numel(Xname)
    %t_str = Xname{i}(end-23:end-9);
    %t(i) = datetime(datenum(t_str,'yyyymmdd_HHMMSS'),'ConvertFrom','datenum');
    load(fullfile(datadir,Xname{i}));
    glint_idx(i) = compute_glint_idx(roi.data_bs,roi.data,3,0);
    if glint_idx(i) > 0.5
        %illu_name = fullfile(savedir,strcat('BR',num2str(i),'.png'));
        %compute_glint_idx(roi.data_bs,roi.data,3,1,illu_name);
        disp(Xname(i));
        %disp(i);
        k=k+1;
    end
end
%t = t';   
glint_idx = glint_idx';


figure;
title(sprintf('N = %u',numel(idx)));
histogram(round(idx),[0.025:0.05:1.025]);
%histogram(compact(idx),[0.025:0.05:1.025]);
%set(gca,'xlim',[0 1]);



%figure; 
%scatter(round,compact);

%% study one particle
%filename = 'part134__CPI__20151112_192003_033_.mat';
%load(fullfile(datadir,filename));
% %%
% M = roi.data_bs;
% M = M(:,:,3); M = 1-M;
% M = double(M);
% Mbin = imfill(roi.data,'holes');
% M(Mbin == 0) = NaN;
% [h,w] = size(M);
% c = [h/2, w/2];
% y = 1:h;
% x = 1:w;
% D = sqrt((y.' - c(1)) .^ 2 + (x - c(2)) .^ 2); 
% 
% dist = D(:);
% intens = M(:);
% intens(intens == 0) = NaN;
% 
% % compute intensity at the center
% dist_thresh = 2.5; % pixels
% dist_c = dist(dist < dist_thresh);
% intens_in = quantile(intens(dist < dist_thresh),0.05);
% 
% % compute average intensity within the "black" part of the particle
% D_min = 0.25 * max(dist);
% D_max = 0.75 * max(dist);
% intens_out = nanmedian(intens(dist >= D_min & dist <= D_max));
% 
% INDEX = (intens_in - intens_out)/(intens_in + intens_out);
% 
% figure; hold on;
% subplot(121); hold on; title(sprintf('IDX = %2.2f',INDEX));
% imshow(roi.data_bs);
% subplot(122);
% plot(dist,intens,'ko');
% xlabel('Distance from centroid');
% ylabel('Intensity [0-255]');


