close all
clear all

load('major_online_Genie_LMMSE_CDLA_varysinr_inr_3.mat','sinr_w','nmse_opt');
nmse_LMMSE=nmse_opt;
load('major_online_baseline_omp_CDLA_varysinr_inr_3.mat','sinr_w','nmse_opt');
omp=nmse_opt;
load('major_online_baseline_vamp_CDLA_varysinr_inr_3.mat','sinr_w','nmse_opt');
vamp=nmse_opt;
load('online_baseline_lasso_CDLA_varysinr_inr_3_v2.mat','sinr_w','nmse_opt');
lasso=nmse_opt;
load('CNN_major_snrvary_online_ici_test_cdlA_varysinr_inr_3.mat','mse_vbioff','mse_nnoff','mse_vbigenie','mse_nngenie');
vbi_off=mse_vbioff;
vbi_genie=mse_vbigenie;
dnn_off=mse_nnoff;
dnn_genie=mse_nngenie;
load('major_ici_online_60k_novbi_sinrvary_3_inr_3_test.mat','mse_nn');
no_vbi=mse_nn;
%%
figure
semilogy(sinr_w,omp,'bv--',sinr_w,lasso,'m^-',sinr_w,vamp,'bx-',sinr_w,no_vbi,'bo-',sinr_w,vbi_off,'m+--',sinr_w,dnn_off,'r^--',...
    sinr_w,vbi_genie,'mx--',sinr_w,nmse_LMMSE,'ko-',sinr_w,dnn_genie,'rv--','LineWidth',1.1);
grid on
hold on
xlabel('SINR (dB)');ylabel('NMSE');
legend('OMP','Modified LASSO','VAMP','S-DSAE (offline)','S-VBI (offline)','Proposed (offline)','S-VBI (online)',...
    'GD-LMMSE','Proposed (online)','FontSize',9,'location','northeast');
%%
% % ylim([6e-2,1e1]);
% text(-5,0.1,'Proposed (offline)','Color','r','FontSize',11);
% text(-5,0.1,'Proposed (online)','Color','r','FontSize',11);
% text(-5,0.1,'S-DSAE (offline)','Color','b','FontSize',11);
% text(-5,0.1,'S-VBI (offline)','Color','m','FontSize',11);
% text(-5,0.1,'S-VBI (online)','Color','m','FontSize',11);
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
% x = [0.1,0.2];y = [0.3,0.4];
% annotation('arrow',x,y,'Color','m','Linewidth',1);
% x = [0.1,0.2];y = [0.3,0.4];
% annotation('arrow',x,y,'Color','m','Linewidth',1);