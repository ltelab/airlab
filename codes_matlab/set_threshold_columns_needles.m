% small script to find a nice threshold between columns and needles
clearvars; close all;

columns_dir = '../training_set/CPI_columns_vs_needles/columns';
needles_dir = '../training_set/CPI_columns_vs_needles/needles';

columns_list = dir(fullfile(columns_dir,'*.mat'));
columns_list = {columns_list.name}';

needles_list = dir(fullfile(needles_dir,'*.mat'));
needles_list = {needles_list.name}';

columns = zeros(numel(columns_list),4);
needles = zeros(numel(needles_list),4);

for i=1:numel(columns_list)
    
    load(fullfile(columns_dir,columns_list{i}));
    
    columns(i,1) = roi.Rect.aspect_ratio; % rect AR
    columns(i,2) = roi.Rect.eccentricity; % rect ecc.
    columns(i,3) = roi.D90.AR;            % D90 AR
    columns(i,4) = roi.E.b/roi.E.a;       % ellipse AR
   
end

for i=1:numel(needles_list)
    
    load(fullfile(needles_dir,needles_list{i}));
    
    needles(i,1) = roi.Rect.aspect_ratio; % rect AR
    needles(i,2) = roi.Rect.eccentricity; % rect ecc.
    needles(i,3) = roi.D90.AR;            % D90 AR
    needles(i,4) = roi.E.b/roi.E.a;       % ellipse AR
    
    
end

c1 = [255,247,188]./255;
c2 = [254,227,145]./255;
c3 = [254,196,79]./255;
c4 = [254,153,41]./255;

hc1 = c4;
hc2 = c2;
halpha = 0.9;

hbins = [0:0.04:1];

figure;
subplot(221); hold on; box on; grid on;
h11 = histogram(columns(:,1),hbins,'Facecolor',hc1,'Facealpha',halpha);
h12 = histogram(needles(:,1),hbins,'Facecolor',hc2,'Facealpha',halpha);
limlim = ylim;
plot([0.26 0.26],[limlim(1) limlim(2)],'r--');
title('Rect. AR');
ylabel('count');

subplot(222); hold on; box on; grid on;
h21 = histogram(columns(:,2),hbins,'Facecolor',hc1,'Facealpha',halpha);
h22 = histogram(needles(:,2),hbins,'Facecolor',hc2,'Facealpha',halpha);
title('Rect. ecc');
ylabel('count');

subplot(223); hold on; box on; grid on;
h31 = histogram(columns(:,3),hbins,'Facecolor',hc1,'Facealpha',halpha);
h32 = histogram(needles(:,3),hbins,'Facecolor',hc2,'Facealpha',halpha);
title('D90 AR');
ylabel('count');

subplot(224); hold on; box on; grid on;
h41 = histogram(columns(:,4),hbins,'Facecolor',hc1,'Facealpha',halpha);
h42 = histogram(needles(:,4),hbins,'Facecolor',hc2,'Facealpha',halpha);
title('Ellipse fit AR');
ylabel('count');



