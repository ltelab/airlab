% plot forward feature selection output
clearvars; close all;

dir1 = 'feat_opt/BER_based/HVPS_1426s_10it_4k.mat';
dir2 = 'feat_opt/BER_based/2DS_4217s_10it_4k.mat';
dir3 = 'feat_opt/BER_based/CPI_2964s_10it_4k.mat';

load(dir1);
featvec_1 = feat_vec_current;
BER_Te_1 = mean_BER_Te;
BER_Tr_1 = mean_BER_Tr;
HSS_Te_1 = mean_kappa_Te;
HSS_Tr_1 = mean_kappa_Tr;

load(dir2);
featvec_2 = feat_vec_current;
BER_Te_2 = mean_BER_Te;
BER_Tr_2 = mean_BER_Tr;
HSS_Te_2 = mean_kappa_Te;
HSS_Tr_2 = mean_kappa_Tr;

load(dir3);
featvec_3 = feat_vec_current;
BER_Te_3 = mean_BER_Te;
BER_Tr_3 = mean_BER_Tr;
HSS_Te_3 = mean_kappa_Te;
HSS_Tr_3 = mean_kappa_Tr;

fv_mat = [sort(featvec_1(1:15)) sort(featvec_2(1:15)) sort(featvec_3(1:15))];



%%


figure; hold on; grid on; box on;
plot(BER_Te_1,'b-','linewidth',1.2);
plot(BER_Tr_1,'b--','linewidth',1.2);
plot(BER_Te_2,'g-','linewidth',1.2);
plot(BER_Tr_2,'g--','linewidth',1.2);
plot(BER_Te_3,'r-','linewidth',1.2);
plot(BER_Tr_3,'r--','linewidth',1.2);
axis([0 30 0 0.4]);
xlabel('# of features used');
ylabel('BER');
legend('HVPS valid','HVPS train','2D-S valid','2D-S train','CPI valid','CPI train');
set(gca,'Fontsize',14);


figure; hold on; grid on; box on;
plot(HSS_Te_1,'b-','linewidth',1.2);
plot(HSS_Tr_1,'b--','linewidth',1.2);
plot(HSS_Te_2,'g-','linewidth',1.2);
plot(HSS_Tr_2,'g--','linewidth',1.2);
plot(HSS_Te_3,'r-','linewidth',1.2);
plot(HSS_Tr_3,'r--','linewidth',1.2);
axis([0 30 0.6 1]);
xlabel('# of features used');
ylabel('HSS');
legend('HVPS valid','HVPS train','2D-S valid','2D-S train','CPI valid','CPI train');
set(gca,'Fontsize',14);