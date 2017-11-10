% small script to analyze cropped images and determine the optimal
% threshold to classify images as truncated
clearvars; close all;

dir_data = '/home/praz/Documents/airlab/training_set/CPI4';
t_str_start = '20100101000000';
t_str_stop = '20200101000000';


[X,Xlab,Xname,Xt] = load_processed_2DS_data(dir_data,t_str_start,t_str_stop,[67 98]);
y = load_2DS_labels(dir_data,t_str_start,t_str_stop);

% remove NaN-values in X,y
idx = find(isnan(sum(X,2)));
X(idx,:) = [];
y(idx) = [];
labels = {'Agg','Col','Gra','Ros','Sph','Pla'};

%%

edge = 0.325;

if 0
    
    figure;
    for i=1:numel(unique(y))

        idx = find(y==i);
        subplot(3,2,i); hold on; grid on; box on;
        histogram(X(idx,1),[0:0.05:1],'Normalization','probability');
        plot([edge edge],[0 1],'r--');
        title(labels{i});
        xlabel('p_{edge}/p_{tot}');

    end

    figure; hold on; box on; grid on;
    idx1 = find(y~=6);
    idx2 = find(y==6);
    histogram(X(idx1,1),[0:0.05:1],'Normalization','probability');
    histogram(X(idx2,1),[0:0.05:1],'Normalization','probability');
    xlabel('p_{edge}/p_{tot}');
    legend('all habits','truncated');
    plot([edge edge],[0 0.5],'r--');

end

figure;
subplot(121);
histogram(X(:,1));
xlabel('p_{edge}/p_{tot}');
subplot(122);
histogram(X(:,2));
xlabel('p_{edge}/frametot');



%histogram(X(:,
