% small piece of code to illustrate the concept of invariant angle and line
% segment (IC-PCA) on different shapes

%clearvars; 
close all;

img_name = 'column.png';
img = imread(img_name); % grayscale or truecolor
data = false(size(img,1),size(img,2));
data(img==0) = true;

P = bwboundaries(data);
P = P{1};
gamma = linspace(0,360,numel(P(:,1)));
Dmax = compute_Dmax(P(:,1),P(:,2),0);


shift_angle = [30, 90, 180];
line_seg_30 = zeros(1,numel(gamma));
line_seg_90 = zeros(1,numel(gamma));
line_seg_180 = zeros(1,numel(gamma));

for i=1:numel(gamma)
    
    gamma_current = gamma(i);
    
    % 30 degrees shift
    tmp_gamma_30 = abs((mod(gamma_current+shift_angle(1),360))-gamma);
    [~,idx_30] = min(tmp_gamma_30);
    A = P(i,:);
    B_30 = P(idx_30,:);
    line_seg_30(i) = sqrt((A(1)-B_30(1))^2 + (A(2)-B_30(2))^2);    
    
    % 90 degrees shift
    tmp_gamma_90 = abs((mod(gamma_current+shift_angle(2),360))-gamma);
    [~,idx_90] = min(tmp_gamma_90);
    A = P(i,:);
    B_90 = P(idx_90,:);
    line_seg_90(i) = sqrt((A(1)-B_90(1))^2 + (A(2)-B_90(2))^2);
    
    % 180 degrees shift
    tmp_gamma_180 = abs((mod(gamma_current+shift_angle(3),360))-gamma);
    [~,idx_180] = min(tmp_gamma_180);
    A = P(i,:);
    B_180 = P(idx_180,:);
    line_seg_180(i) = sqrt((A(1)-B_180(1))^2 + (A(2)-B_180(2))^2);
    
end
    
line_seg_30N = line_seg_30 ./ Dmax;
line_seg_90N = line_seg_90 ./ Dmax;
line_seg_180N = line_seg_180 ./ Dmax;


% compute ICPCA stuff
outputDir = '/ltedata/MASC/OAP/airlab/codes_matlab/invariant_angle';
       % icpca.plots = false;
       % icpca.debug = false;
       % icpca.enableFilters = false;
       % if process.icpca_default % 1:use raw image as input | 0:use DS.data as input 
[icpca.x,icpca.y,icpca.amax,icpca.area_no_holes,icpca.name,icpca.cout,icpca.centroid,icpca.adiff,icpca.imsize,icpca.edgepx,icpca.xtremepx] = get_perimeter(img_name, outputDir, true, false, false);
 %       else
  %          [icpca.x,icpca.y,icpca.amax,icpca.area_no_holes,icpca.name,icpca.cout,icpca.centroid,icpca.adiff,icpca.imsize,icpca.edgepx,icpca.xtremepx] = get_perimeter(icpca.filePath, icpca.outputDir, icpca.plots, icpca.debug, icpca.enableFilters,double(DS.data));
  %      end
bullet = get_parameters(img_name, outputDir, outputDir, outputDir, 1, 1, icpca.x, icpca.y, icpca.amax, icpca.area_no_holes, icpca.name, icpca.cout, icpca.centroid, icpca.adiff, icpca.imsize, icpca.edgepx, icpca.xtremepx, 0, outputDir);
%disp(DS.icpca.b.jana_ka);
%disp(DS.icpca.b.jana_kov);



figure(1);
imshow(~data);hold on;
set(gca,'Ydir','reverse');
hold on;
plot(P(:,2),P(:,1),'r.-','linewidth',2);
plot(P(1,2),P(1,1),'o','MarkerFaceColor','w','MarkerEdgeColor','k','MarkerSize',10,'linewidth',2);
idx = find(gamma > 30,1,'first');
plot(P(idx,2),P(idx,1),'o','MarkerFaceColor',[222 45 38]./255,'MarkerEdgeColor','k','MarkerSize',7,'linewidth',2);
idx = find(gamma > 90,1,'first');
plot(P(idx,2),P(idx,1),'o','MarkerFaceColor',[65,171,93]./255,'MarkerEdgeColor','k','MarkerSize',7,'linewidth',2);
idx = find(gamma > 180,1,'first');
plot(P(idx,2),P(idx,1),'o','MarkerFaceColor',[37,52,148]./255,'MarkerEdgeColor','k','MarkerSize',7,'linewidth',2);
figure;
hold on; box on;
plot(gamma,line_seg_180N,'-','color',[37,52,148]./255,'linewidth',1.5);
plot(gamma,line_seg_90N,'-','color',[65,171,93]./255,'linewidth',1.5);
plot(gamma,line_seg_30N,'-','color',[222 45 38]./255,'linewidth',1.5);
set(gca,'ylim',[0 1.1]);
set(gca,'xlim',[0 360]);
set(gca,'XTick',[0 45 90 135 180 225 270 315 360]);
set(gca,'YTick',[0 0.5 1]);
xlabel('shift from P_0');
ylabel('s(\gamma)/D_{max}');
set(gca,'Fontsize',12);

%% 
% clearvars -except circle column bullet
% 
% B = [circle.b.jana_kov; column.b.jana_kov; bullet.b.jana_kov];
% figure;
% bar(B');


% 
% cov(gamma(90),30); %30 = + 27 shift in gamme
% idx_shift30 = find(gamma >= 30,1,'first');
% [autocov,lag] = xcov(line_seg_180N);
% autocov = autocov./max(autocov);
% figure(9);
% plot(lag,autocov); hold on;
% plot([lag(numel(line_seg_180N)+idx_shift30) lag(numel(line_seg_180N)+idx_shift30)],[0 1],'r-'); hold off;

