% small code to illustrate best MLR hyperparameters combinations
clearvars; close all;
load('./grid_search/CPI_MLR_N2964.mat');

% BER
tmp_min = min(BER_Te(:));
[idx1,idx2] = find(BER_Te == tmp_min);
step_min_BER = step_vec(idx1);
lambda_min_BER = lambda_vec(idx2);

figure; hold on;
contourf(log10(step_vec),log10(lambda_vec),100.*BER_Te',1000,'Edgecolor','None');
h1 = plot3(log10(step_min_BER),log10(lambda_min_BER),100,'ro','MarkerFaceColor','r','MarkerEdgeColor','k');
l1 = 'MIN';
c = colorbar;
c.Label.String = 'BER [%]';
box on;
xlabel('log_{10}(gradient step)');
ylabel('log_{10}(\lambda)');
legend(h1,l1);
title('BER');

% HSS
tmp_min = max(kappa_Te(:));
[idx1,idx2] = find(kappa_Te == tmp_min);
step_max_kappa = step_vec(idx1);
lambda_max_kappa = lambda_vec(idx2);

figure; hold on;
contourf(log10(step_vec),log10(lambda_vec),100.*kappa_Te',1000,'Edgecolor','None');
h1 = plot3(log10(step_max_kappa),log10(lambda_max_kappa),100,'ro','MarkerFaceColor','r','MarkerEdgeColor','k');
l1 = 'MAX';
c = colorbar;
c.Label.String = 'kappa [%]';
box on;
xlabel('log_{10}(gradient step)');
ylabel('log_{10}(\lambda)');
legend(h1,l1);
title('Cohen''s kappa');

% OA
tmp_min = max(OA_Te(:));
[idx1,idx2] = find(OA_Te == tmp_min);
step_max_OA = step_vec(idx1);
lambda_max_OA = lambda_vec(idx2);

figure; hold on;
contourf(log10(step_vec),log10(lambda_vec),100.*OA_Te',1000,'Edgecolor','None');
h1 = plot3(log10(step_max_OA),log10(lambda_max_OA),linspace(100,100,numel(step_max_OA)),'ro','MarkerFaceColor','r','MarkerEdgeColor','k');
l1 = 'MAX';
c = colorbar;
c.Label.String = 'OA [%]';
box on;
xlabel('log_{10}(gradient step)');
ylabel('log_{10}(\lambda)');
legend(h1,l1);
title('Overall Accuracy');



% subplot(132);
% contourf(log10(step_vec),log10(lambda_vec),kappa_Te',1000,'Edgecolor','None');
% colorbar;
% subplot(133);
% contourf(log10(step_vec),log10(lambda_vec),OA_Te',1000,'Edgecolor','None');
% colorbar;


