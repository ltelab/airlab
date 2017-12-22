% investigate the value of the new IC-PCA descriptors on the classification
clear all; close all;

dir_data = '/home/praz/Documents/Cloud_Probes/training_set/2DS';
t_str_start = '20150101000000';
t_str_stop  = '20180101000000';

[X,Xlab,Xname,Xt] = load_processed_2DS_data(dir_data,t_str_start,t_str_stop);
y = load_2DS_labels(dir_data,t_str_start,t_str_stop);

%% do stuff

X_ic = X(:,70:end);
labels = {'Agg','Col','Gra','Ros','Sph','Oth'};
Xlab_ic = Xlab(70:end);


for i=1:2:size(X_ic,2)-1

    figure;
    scatterhist(X_ic(:,i),X_ic(:,i+1),'Group',y);
    xlabel(sprintf('feat. #%s',Xlab_ic{i}));
    ylabel(sprintf('feat. #%s',Xlab_ic{i+1}));
    
end

figure;
scatterhist(X_ic(:,6),X_ic(:,8),'Group',y,'Style','bar','MarkerSize',3);
xlabel(sprintf('feat. #%s',Xlab_ic{6}));
ylabel(sprintf('feat. #%s',Xlab_ic{8}));