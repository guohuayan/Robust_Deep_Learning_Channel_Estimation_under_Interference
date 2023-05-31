close all
clear all

addpath('./matlab_lib');
load('major_online_vbitrain_CDLA_sub60_test_1e5.mat','cov_H_mix','power_d','Fm','Fa','Fi','R_w');
%%
sinr_w = linspace(-10,0,5);
nmse_opt=zeros(1,length(sinr_w));
%%
for s0=1:length(sinr_w)
    %%
    if s0==1
        load('online_test_CDLA_60k_sinr_n10_Inr_3.mat','HLS_delay_w','Ct_int_w','noise_p','Ht_freq_w','Cu_inv_w');
    elseif s0==2
        load('online_test_CDLA_60k_sinr_n7.5_Inr_3.mat','HLS_delay_w','Ct_int_w','noise_p','Ht_freq_w','Cu_inv_w');
    elseif s0==3
        load('online_test_CDLA_60k_sinr_n5_Inr_3.mat','HLS_delay_w','Ct_int_w','noise_p','Ht_freq_w','Cu_inv_w');
    elseif s0==4
        load('online_test_CDLA_60k_sinr_n2.5_Inr_3.mat','HLS_delay_w','Ct_int_w','noise_p','Ht_freq_w','Cu_inv_w');
    elseif s0==5
        load('online_test_CDLA_60k_sinr_n0_Inr_3.mat','HLS_delay_w','Ct_int_w','noise_p','Ht_freq_w','Cu_inv_w');
    end
    %%
    NSlots=1000;
    HLS_delay_w=HLS_delay_w(:,:,1:NSlots);
    Ct_int_w=Ct_int_w(:,:,1:NSlots);
    Ht_freq_w=Ht_freq_w(:,:,1:NSlots);
    %%
    Nd=size(Ht_freq_w,1);
    M=size(Ht_freq_w,2);
    %%
    H_delay_w=HLS_delay_w;
    n_delay=size(H_delay_w,1);
    H_delay_obs=permute(H_delay_w,[2,1,3]); %% M,n_delay,NSlots
    %%       
    Hopt_delay_w=zeros(M,n_delay,NSlots);
    parfor tau0=1:NSlots
        Ct_int=Ct_int_w(:,:,tau0);
        H_obs=H_delay_obs(:,:,tau0);
        %%
        if mod(tau0,100)==0
            fprintf('t0=%d\n',tau0);
        end
        %%
%         lambda=sqrt(noise_p);
        Cu=Ct_int+noise_p*eye(M);
        H_e = cvx_lasso_ici(H_obs,Cu);
        %%
        Hopt_delay_w(:,:,tau0)=H_e;
    end
    H_delay_w=permute(Hopt_delay_w,[2,1,3]); %% n_delay,M,NSlots
    H_est_w=Fm*reshape(H_delay_w,[n_delay,M*NSlots]);
    H_est_w=reshape(H_est_w,[Nd,M,NSlots]);
    err=Ht_freq_w-H_est_w;
    NMSE=sum(abs(err).^2,[1,2])./sum(abs(Ht_freq_w).^2,[1,2]);
    nmse_opt(s0)=mean(NMSE(:));
    mean(NMSE(:))
end
figure
semilogy(sinr_w,nmse_opt,'b+-');
grid on
%%
save('online_baseline_lasso_CDLA_varysinr_inr_3_v2.mat','sinr_w','nmse_opt');