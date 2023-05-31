function [power_d_o,cov_H_mix_o,Hoe_delay_w] = online_VBI_training_func(speedin,He_delay_buff,Cd_int_buff,power_d_o,cov_H_mix_o,noise_d,io)
    [n_delay,M,num_pilot,Lwin]=size(He_delay_buff);
    num_train=num_pilot*Lwin;
    Ht_delay_o=reshape(He_delay_buff,[n_delay,M,num_train]); %% n_delay,M,num_train
    Cwt_int_o=reshape(Cd_int_buff,[M,M,num_train]);
    %%
    Ho_delay_obs=permute(Ht_delay_o,[2,3,1]); %% M,num_train,n_delay
    Cov_Hef=zeros(M,M,n_delay);
    for t0=1:n_delay
        Cov_Hef(:,:,t0)=power_d_o(t0)*cov_H_mix_o;
    end
    %%
    Hoe_delay_w=zeros(M,num_train,n_delay);
    Cov_x=zeros(M,M,n_delay);
    for j0=1:num_train
        Co_int=Cwt_int_o(:,:,j0)+eye(M)*noise_d;
        Cov_xt=zeros(M,M,n_delay);
        for t0=1:n_delay
            Ch=Cov_Hef(:,:,t0);
            filter=Ch/(Ch+Co_int);
            Ht_e=Ho_delay_obs(:,j0,t0);
            Ht_e=filter*Ht_e;
            Hoe_delay_w(:,j0,t0)=Ht_e;
            Cx=inv(inv(Ch)+inv(Co_int))+(Ht_e*Ht_e');
            Cov_xt(:,:,t0)=Cx;
        end
        Cov_x=Cov_x+Cov_xt;
    end
    %%
%     Hoe_delay_w=permute(Hoe_delay_w,[2,3,1]); %% Lwin,n_delay,M
    Hoe_delay_w=permute(Hoe_delay_w,[3,1,2]); %% n_delay,M,num_train
    Hoe_delay_w=reshape(Hoe_delay_w,[n_delay,M,num_pilot,Lwin]);
    Hoe_delay_w=Hoe_delay_w(:,:,:,end);
    Hoe_delay_w=permute(Hoe_delay_w,[1,3,2]); %% n_delay,num_pilot,M
    %%
    Cov_x=Cov_x./num_train;
    npower_d=zeros(n_delay,1);
%         for t0=1:n_delay
%             Cx=Cov_x(:,:,t0);
%             npower_d(t0)=trace(Cx)/M;
%         end
    Covi=inv(cov_H_mix_o);
    for t0=1:n_delay
        Cx=Cov_x(:,:,t0);
        npower_d(t0)=trace(Cx*Covi)/M;
    end
    ncov_H_mix=zeros(M,M);
    for t0=1:n_delay
        Cx=Cov_x(:,:,t0)./npower_d(t0);
        ncov_H_mix=ncov_H_mix+Cx;
    end
    ncov_H_mix=ncov_H_mix./trace(ncov_H_mix)*M;
    %%
%         dpower_d=npower_d-power_d_o;
    dpower_d=1./npower_d-1./power_d_o;
    dcov_H_mix=inv(ncov_H_mix)-inv(cov_H_mix_o);
    %%
    if speedin == 3
        rho=10/(200+io); %% speed=3
    else
        rho=2/(100+io)^0.6; %% speed=30 60 120
    end
    power_d_o=1./power_d_o+rho*dpower_d;
    power_d_o=1./power_d_o;
    cov_H_mix_o=inv(cov_H_mix_o)+rho*dcov_H_mix;
    cov_H_mix_o=inv(cov_H_mix_o);
end

