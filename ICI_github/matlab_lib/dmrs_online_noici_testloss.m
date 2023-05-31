function nmse = dmrs_online_noici_testloss(power_d,cov_H_mix,Ht_delay_w,noise_p,Fm,Ht_freq_w)
    Ht_delay_obs=permute(Ht_delay_w,[2,3,1]); %% M,NSlots,n_delay
    [M,NSlots,n_delay]=size(Ht_delay_obs);
    Cov_Hef=zeros(M,M,n_delay);
    for t0=1:n_delay
        Cov_Hef(:,:,t0)=power_d(t0)*cov_H_mix;
    end
    %%       
    Hte_delay_w=zeros(M,NSlots,n_delay);
%     est_var=zeros(n_delay,1);
    for t0=1:n_delay
        Ch=Cov_Hef(:,:,t0);
        filter=Ch/(Ch+eye(M)*noise_p);
        Ht_e=Ht_delay_obs(:,:,t0);
        Ht_e=filter*Ht_e;
        Hte_delay_w(:,:,t0)=Ht_e;
    end
    Ht_delay_w=permute(Hte_delay_w,[3,1,2]); %% n_delay,M,NSlots
    NMSE_vector=zeros(NSlots,1);
    for i0=1:NSlots
        Ht_delay=Ht_delay_w(:,:,i0);
        Ht_est=Fm*Ht_delay;
        Ht_freq=Ht_freq_w(:,:,i0);
        err=Ht_freq-Ht_est;
        NMSE_vector(i0)=sum(abs(err(:)).^2)/sum(abs(Ht_freq(:)).^2);
    end
    nmse=mean(NMSE_vector(:));
end

