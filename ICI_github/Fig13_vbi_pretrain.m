close all
clear all

addpath('./matlab_lib');
load('major_off_test_mixmodel_sub60_1e3.mat','H_home_batch');
load('offline_vbitrain_major_sub60_1e5.mat','cov_H_mix','power_d','Fm','Fa','Fi','R_w');
%%
% load('dmrs_channel_256_home_mixmodel_1e3.mat','H_home_batch');
% load('offline_256_vbitrain_1e5.mat','cov_H_mix','power_d','Fm','Fa','Fi','R_w');
%%
[Nd,M,NSlots]=size(H_home_batch);
nfft=Nd;
%%
F_pilot=generate_ZC_pilot(Nd,1,0);
%%
snr_w=linspace(-10,0,5);
% snr_w=10;
nmse_vbi=zeros(1,length(snr_w));
nmse_opt=zeros(1,length(snr_w));
%%
for s0=1:length(snr_w)
    SNRdB=snr_w(s0);
    fprintf('snr=%d\n',SNRdB);
    Pt = 10^(SNRdB/20);
    H_freq_w=zeros(Nd,M,NSlots);
    H_LS_w=zeros(Nd,M,NSlots);
    for i0=1:NSlots
        H_freq=H_home_batch(:,:,i0);
        Rx=diag(F_pilot)*H_freq+AWGN(size(H_freq))./Pt;
%         Rx=diag(F_pilot)*H_freq;
        H_est=diag(F_pilot')*Rx;
%         H_est=H_freq+AWGN(size(H_freq))./Pt;
        %%
        H_freq_w(:,:,i0)=H_freq;
        H_LS_w(:,:,i0)=H_est;
    end 
    %%
    noise_p=(1/Pt)^2;
    n_delay=size(Fm,2);
    H_delay_w=zeros(n_delay,M,NSlots);
    for i0=1:NSlots
        %%
        H_est=H_LS_w(:,:,i0);
        H_delay=Fm'*H_est;
        H_delay_w(:,:,i0)=H_delay;
    end
    H_delay_obs=permute(H_delay_w,[2,3,1]); %% M,NSlots,n_delay
    Cov_Hef=zeros(M,M,n_delay);
    for t0=1:n_delay
        Cov_Hef(:,:,t0)=power_d(t0)*cov_H_mix;
    end
    %%       
    He_delay_w=zeros(M,NSlots,n_delay);
    Hopt_delay_w=zeros(M,NSlots,n_delay);
    for t0=1:n_delay
        H_obs=H_delay_obs(:,:,t0);
        Ch=Cov_Hef(:,:,t0);
        filter=Ch/(Ch+eye(M)*noise_p);
        H_e=filter*H_obs;
        He_delay_w(:,:,t0)=H_e;
        Ch=R_w(:,:,t0);
        filter=Ch/(Ch+eye(M)*noise_p);
        H_e=filter*H_obs;
        Hopt_delay_w(:,:,t0)=H_e;
    end
    H_delay_w=permute(He_delay_w,[3,1,2]); %% n_delay,M,NSlots
    H_est_w=Fm*reshape(H_delay_w,[n_delay,M*NSlots]);
    H_est_w=reshape(H_est_w,[Nd,M,NSlots]);
    err=H_freq_w-H_est_w;
    NMSE=sum(abs(err).^2,[1,2])./sum(abs(H_freq_w).^2,[1,2]);
    nmse_vbi(s0)=mean(NMSE(:));
    H_delay_w=permute(Hopt_delay_w,[3,1,2]); %% n_delay,M,NSlots
    H_est_w=Fm*reshape(H_delay_w,[n_delay,M*NSlots]);
    H_est_w=reshape(H_est_w,[Nd,M,NSlots]);
    err=H_freq_w-H_est_w;
    NMSE=sum(abs(err).^2,[1,2])./sum(abs(H_freq_w).^2,[1,2]);
    nmse_opt(s0)=mean(NMSE(:));
end
figure
semilogy(snr_w,nmse_vbi,'ro-',snr_w,nmse_opt,'b+-');
grid on
save('major_ici_free_60k_vbi_pretrain_n10_0.mat','snr_w','nmse_opt','nmse_vbi');