% % script to compute classification PSD
%close all; 
close all; clearvars;

reload_data = 0;
load('flight2_00h47.mat');

flight_name = ' flight 2 20151201 00h47';
probe_res = [10 150 2.3]; % in um/pixel
probe_name = {'2DS','HVPS','CPI'};
 

%% inputs
datadir{1} = '/ltedata/MASC/OAP/OAP_flight_data/20151201_WF/2DS/00h32_proc'; % 2DS
datadir{2} = '/ltedata/MASC/OAP/OAP_flight_data/20151201_WF/HVPS/00h32_proc'; % HVPS 
datadir{3} = '/ltedata/MASC/OAP/OAP_flight_data/20151201_WF/CPI/00h30-00h35_proc'; % CPI

if reload_data
    
    % initilalization
        Nini = zeros(3,1);
        Nkept = zeros(3,1);
        Dmax = cell(3,1);
        area = cell(3,1);
        labelID = cell(3,1);
        MatP = cell(3,1);
    

    for i=1:numel(datadir)

        [X, Xlab, Xname, Xt, Xfullprob] = load_processed_2DS_labels(datadir{i},[],[],true);
        y = X(:,1);%load_2DS_labels(datadir{i});
        Nini(i) = numel(y);
        Dmax_tmp = X(:,3);
        area_tmp = X(:,2);
        % Dmax estimated from area for particles too small
        Dmax_tmp(isnan(Dmax_tmp) & ~isnan(area_tmp)) = sqrt(area_tmp(isnan(Dmax_tmp) & ~isnan(area_tmp)));
        % store into cells
        Dmax{i} = Dmax_tmp;
        area{i} = area_tmp;
        Dmax{i} = Dmax{i} * probe_res(i) / 1000; % conversion to mm
        labelID{i} = y;
        Nkept(i) = sum(labelID{i}>0);
        MatP{i} = Xfullprob;
        
        % filter labelID < 0 (noisy imgs, errors, ...)
        MatP{i}(labelID{i}<0,:) = [];
        area{i}(labelID{i}<0) = [];
        Dmax{i}(labelID{i}<0) = [];
        labelID{i}(labelID{i}<0) = [];

        fprintf('\n*** %s cloud probe data loaded *****\n',probe_name{i});
        fprintf('Number of images found : %u \n',Nini(i));
        fprintf('Number of images kept  : %u \n',Nkept(i));
        fprintf('************************************\n\n');    

    end
    
end


%% Classification scheme
% -2 = noise
% -1 = error
%  1 = Agg
%  2 = Col
%  3 = Compact
%  4 = BulRose
%  5 = Quasi-sphere
%  6 = Planar (or Other for 2D-S)
%  7 = Planar for 2D-S
%  8 = Small
% 10 = Truncated
% close all;
%
% new scheme:
% 1 = Agg
% 2 = Col
% 3 = Compact
% 4 = BulRose
% 5 = QS
% 6 = Planar
% 7 = Small
% 8 = Trunc

% some specs
merge_small_trunc = true;
N_class = 10;
x_bins = (0:100:6000)/1000;
PSD_matrix = cell(3,1);
fig_colors = jet(N_class);
labels = cell(3,1);

% filter too small particles
Dlim = [2.5 3.74165 1] .* probe_res ./ 1000;
for i=1:3
    
    area{i}(Dmax{i} < Dlim(i)) = [];
    labelID{i}(Dmax{i} < Dlim(i)) = [];
    %MatP{i}(Dmax{i} < Dlim(i),:) = [];
    Dmax{i}(Dmax{i} < Dlim(i)) = [];
 
end

% filter low probability
if 0
    Plim = [0.95 0 0.9];
    for i=1:3

        prob_tmp = MatP{i};
        maxprob = max(prob_tmp')';
        N_prefilt = numel(area{i});
        area{i}(maxprob < Plim(i)) = [];
        N_postfilt = numel(area{i});
        labelID{i}(maxprob < Plim(i)) = [];
        Dmax{i}(maxprob < Plim(i)) = [];
        fprintf('Prob filtering removed %u images for %s (%2.1f%%)\n',N_postfilt,probe_name{i},(1-N_postfilt/N_prefilt)*100);

    end
end



for i=1:3

    PSD_matrix{i} = zeros(N_class,numel(x_bins));
    
    for j=1:N_class   
        Dmax_class = Dmax{i}(labelID{i}==j);
        bins_class = hist(Dmax_class,x_bins);
        PSD_matrix{i}(j,:) = bins_class;
    
    end

    % merge classes
    if strcmp(probe_name{i},'CPI')
        
        if merge_small_trunc
            PSD_matrix{i}(7,:) = PSD_matrix{i}(8,:) + PSD_matrix{i}(10,:);
            PSD_matrix{i}(8:end,:) = [];
            labels{i} = {'AG','CC','CP','BR','QS','PC','SP+Trunc'};
        else
            PSD_matrix{i}(7,:) = PSD_matrix{i}(8,:);
            PSD_matrix{i}(8,:) = PSD_matrix{i}(10,:);
            PSD_matrix{i}(9:end,:) = [];
            labels{i} = {'AG','CC','CP','BR','QS','PC','SP','Trunc'};
        end
  
    elseif strcmp(probe_name{i},'HVPS')
        
        if merge_small_trunc
            PSD_matrix{i}(7,:) = PSD_matrix{i}(8,:) + PSD_matrix{i}(10,:);
            PSD_matrix{i}(8:end,:) = [];
            labels{i} = {'AG','CC','CP','BR','QS','PC','SP+Trunc'}; 
        else
            PSD_matrix{i}(7,:) = PSD_matrix{i}(8,:);
            PSD_matrix{i}(8,:) = PSD_matrix{i}(10,:);
            PSD_matrix{i}(9:end,:) = [];
            labels{i} = {'AG','CC','CP','BR','QS','PC','SP','Trunc'}; 
        end
    
    elseif strcmp(probe_name{i},'2DS')
        
        others = PSD_matrix{i}(6,:);
        planar = PSD_matrix{i}(7,:);
        small = PSD_matrix{i}(8,:);
        trunc = PSD_matrix{i}(10,:);
        PSD_matrix{i}(6,:) = planar;
        
        if merge_small_trunc
            PSD_matrix{i}(7,:) = small + others + trunc;
            PSD_matrix{i}(8:end,:) = [];
            labels{i} = {'AG','CC','CP','BR','QS','PC','SP+Trunc'};
        else
            PSD_matrix{i}(7,:) = small;
            PSD_matrix{i}(8,:) = others + trunc;
            PSD_matrix{i}(9:end,:) = [];
            labels{i} = {'AG','CC','CP','BR','QS','PC','SP','Trunc'};
        end       
        
    end
    
    PSD_matrix{i} = PSD_matrix{i}';
        
        
end

%% illustration
if merge_small_trunc
    N_class = 7;
    N_class_pie = 7; % 6 = ignoring small/trunc
else
    N_class = 8;
    N_class_pie = 8;
end

for i=1:3
  
    % just the PSD shape
%     figure; hold on; box on; grid on;
%     histogram(Dmax{i},x_bins);
%     title(sprintf('PSD %s %s',probe_name{i},flight_name));
%     xlabel('[mm]');
%     ylabel('Count');

    % PSD by habit
    figure('units','pixels','Position',[100 100 762 638]); hold on; box on; grid on;
    title(sprintf('PSD %s %s',probe_name{i},flight_name));
    b = bar(x_bins,PSD_matrix{i},'stacked');
    for j=1:numel(N_class)
        b(j).FaceColor = fig_colors(j,:);
    end
    legend(labels{i});
    xlabel('[mm]');
    ylabel('Count');
    set(gca,'xlim',[-0.07 4]);
    
end

% pie charts
pie_count = zeros(3,N_class_pie);
for i=1:3
    totals = sum(PSD_matrix{i});
    for j=1:N_class_pie
        pie_count(i,j) = totals(j);
        pie_perc(i,j) = round(100*totals(j)/sum(totals));
    end
end

% add an arbitrary small numbers to pie_count element == 0 (for display)
pie_count(pie_count==0) = 0.01;

if N_class_pie == 7
    pie_labels = {'AG','CC','CP','BR','QS','PC','SP+Trunc'};
else
    pie_labels = {'AG','CC','CP','BR','QS','PC','SP','Trunc'};    
end

slices_labels = {};
figure('units','pixels','Position',[0 0 2000 638]);
for i=1:3
       
    for j=1:N_class_pie
        slices_labels{j} = strcat(pie_labels{j},' ',num2str(pie_perc(i,j)),'%');
    end
    
    subplot(1,3,i);
    %slices_labels = strcat(pie_labels,pie_spaces,num2str(pie_perc(i,:)),pie_symbols);
    pie(pie_count(i,:),slices_labels);
end




