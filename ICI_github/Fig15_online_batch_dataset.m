close all
clear all

% load('test_CDLA_3e3.mat','H_home_batch','H_int_batch');
% load('online_vbitrain_CDLA_test_1e5.mat','cov_H_mix','power_d','Fm','Fa','Fi','R_w');
load('major_test_CDLA_sub60_3e3.mat','H_home_batch','H_int_batch');
load('major_online_vbitrain_CDLA_sub60_test_1e5.mat','cov_H_mix','power_d','Fm','Fa','Fi','R_w');
%%
H_home_test=H_home_batch;
H_int_test=H_int_batch;
%%
%%
[Nd,M,NtSlots]=size(H_home_test);
%%
F_pilot=generate_ZC_pilot(Nd,1,0);
%% interference cell
Ci=size(H_int_test,3);
Intercell_pilot=cell(Ci,1);
cv=physconst('LightSpeed');
for c0=1:Ci
    %% pilot
    F_pilot_int_base=generate_ZC_pilot(Nd,c0+2,0);
    Intercell_pilot{c0}.F_pilot_int=F_pilot_int_base;
end
%%
% snr_w=[-15, -10, -5, 0, 5];
sinr_w = linspace(-10,0,5); % -10, -7.5 -5 -2.5 0
Inr= 3;
Pin=10^(Inr/10);
len_snr=length(sinr_w);
% channel_sampe=cell(len_snr,1);
%%
for s0=1:len_snr
    SNRdB=sinr_w(s0);
    %%
    fprintf('snr=%d\n',SNRdB);
    %%
    P_ipn=10^(-SNRdB/10);
    Pn = P_ipn./(1+Pin);
    Pint=Pn*Pin;
    %%
    Pt=1./sqrt(Pn);
    Pint=sqrt(Pint);
    %%
%     Pt = 10^(SNRdB/20);
%     Pint=sqrt(10^(Inr/10))./Pt;
    Ht_freq_w=zeros(Nd,M,NtSlots);
    Ht_LS_w=zeros(Nd,M,NtSlots);
    Ct_int_w=zeros(M,M,NtSlots);
    for i0=1:NtSlots
        Ht_freq=H_home_test(:,:,i0);
        Rx=diag(F_pilot)*Ht_freq+AWGN(size(Ht_freq))./Pt;
        for c0=1:Ci
            pilot=Intercell_pilot{c0}.F_pilot_int;
            H_int_u=H_int_test(:,:,c0,i0);
            Rx=Rx+Pint.*diag(pilot)*H_int_u;
        end
        Ht_est=diag(F_pilot')*Rx;
        %%
        He_int=Fi'*Ht_est;
        noise_var=1/Pt^2;
%         C_int = ici_cov_est(He_int,noise_var,nfft,Np,M);
        C_int = ici_cov_est_ofdm(He_int,noise_var,M);
        %%
        Ht_freq_w(:,:,i0)=Ht_freq;
        Ht_LS_w(:,:,i0)=Ht_est;
        Ct_int_w(:,:,i0)=C_int;
    end 
    n_delay=size(Fm,2);
    HLS_delay_w=zeros(n_delay,M,NtSlots);
    for i0=1:NtSlots
        %%
        Ht_est=Ht_LS_w(:,:,i0);
        Ht_delay=Fm'*Ht_est;
        HLS_delay_w(:,:,i0)=Ht_delay;
    end
    %%
    Cu_inv_w=zeros(NtSlots,M,M);
    for i0=1:NtSlots
        Cint=Ct_int_w(:,:,i0)+(1/Pt)^2*eye(M);
        Cu_inv_w(i0,:,:)=inv(Cint);
    end
    %% delay domain variable
    noise_p=(1/Pt)^2;
%     Ct_int_w=Ct_int_w;
    %%
    if s0==1
        save('online_test_CDLA_60k_sinr_n10_Inr_3.mat','HLS_delay_w','Ct_int_w','noise_p','Ht_freq_w','Cu_inv_w');
    elseif s0==2
        save('online_test_CDLA_60k_sinr_n7.5_Inr_3.mat','HLS_delay_w','Ct_int_w','noise_p','Ht_freq_w','Cu_inv_w');
    elseif s0==3
        save('online_test_CDLA_60k_sinr_n5_Inr_3.mat','HLS_delay_w','Ct_int_w','noise_p','Ht_freq_w','Cu_inv_w');
    elseif s0==4
        save('online_test_CDLA_60k_sinr_n2.5_Inr_3.mat','HLS_delay_w','Ct_int_w','noise_p','Ht_freq_w','Cu_inv_w');
    elseif s0==5
        save('online_test_CDLA_60k_sinr_n0_Inr_3.mat','HLS_delay_w','Ct_int_w','noise_p','Ht_freq_w','Cu_inv_w');
    end
end
%%






