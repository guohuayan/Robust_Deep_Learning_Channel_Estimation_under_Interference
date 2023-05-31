close all
clear all

load('major_ici_free_60k_vbi_pretrain_n10_0.mat','snr_w','nmse_opt','nmse_vbi');
vbi_off=nmse_vbi;
LMMSE=nmse_opt;
load('major_ici_free_60k_cnn_offtest.mat', 'snr_w','test_loss');
CNN_prop=test_loss;
load('major_ici_free_60k_cnn_off_novbi.mat', 'snr_w','test_loss');
CNN_novbi=test_loss;
%%
load('major_ici_free_60k_omp_n10_0.mat','snr_w','nmse_vbi');
omp=nmse_vbi;
load('major_ici_free_60k_lasso_n10_0.mat','snr_w','nmse_vbi');
lasso=nmse_vbi;
load('major_ici_free_60k_vamp_n10_0.mat','snr_w','nmse_vbi');
vamp=nmse_vbi;
% close all
figure
semilogy(snr_w,omp,'b^-',snr_w,lasso,'mv-',snr_w,vamp,'bx-',snr_w,vbi_off,'r^--',snr_w,LMMSE,'kx-',...
    snr_w,CNN_novbi,'bo-',snr_w,CNN_prop,'r+--','LineWidth',1.1);
grid on
hold on
xlabel('SNR (dB)');ylabel('NMSE');
legend('OMP','Modified LASSO','VAMP','Proposed (S-VBI)','GD-LMMSE','Proposed (S-DSAE)','Proposed','FontSize',9,'location','southwest');
% legend('OMP','VAMP','Proposed (VBI only)','GD-LMMSE','Proposed (DSAE only)','Proposed','FontSize',10,'location','northeast');
% xlim([-15,5]);
ylim([9e-3,5.1]);
% text(-5,0.1,'Proposed','Color','r','FontSize',11);
% text(-5,0.1,'Proposed (S-DSAE)','Color','b','FontSize',11);
% text(-5,0.1,'Proposed (S-VBI)','Color','r','FontSize',11);
% text(-5,0.1,'GD-LMMSE','Color','k','FontSize',11);
% text(-5,0.8,'OMP','Color','b','FontSize',11);
% text(-5,0.8,'Modified LASSO','Color','m','FontSize',11);
% text(-5,0.8,'VAMP','Color','b','FontSize',11);
% x = [0.3,0.5];y = [0.6,0.5];
% annotation('arrow',x,y,'Color','r','Linewidth',1);
% x = [0.3,0.5];y = [0.3,0.4];
% annotation('arrow',x,y,'Color','b','Linewidth',1);
% x = [0.3,0.5];y = [0.6,0.5];
% annotation('arrow',x,y,'Color','r','Linewidth',1);
% x = [0.3,0.5];y = [0.3,0.4];
% annotation('arrow',x,y,'Color','b','Linewidth',1);
% x = [0.2,0.6];y = [0.3,0.4];
% annotation('arrow',x,y,'Color','b','Linewidth',1);
% x = [0.6,0.7];y = [0.3,0.4];
% annotation('arrow',x,y,'Color','k','Linewidth',1);
% x = [0.1,0.2];y = [0.3,0.4];
% annotation('arrow',x,y,'Color','m','Linewidth',1);
%%