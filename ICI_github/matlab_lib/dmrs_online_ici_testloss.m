function nmse = dmrs_online_ici_testloss(power_d,cov_H_mix,Ht_delay_w,Cwt_int_w,noise_p,Fm,Ht_freq_w)
    H_delay_obs=permute(Ht_delay_w,[2,3,1]); %% M,NSlots,n_delay
    [M,NSlots,n_delay]=size(H_delay_obs);
    Cov_Hef=zeros(M,M,n_delay);
    for t0=1:n_delay
        Cov_Hef(:,:,t0)=power_d(t0)*cov_H_mix;
    end
    He_delay_w=zeros(M,NSlots,n_delay);
    for n0=1:NSlots
        C_int_tmp=Cwt_int_w(:,:,n0);
%         C_int_tmp=0;
        for t0=1:n_delay
            Ch=Cov_Hef(:,:,t0);
            filter=Ch/(Ch+eye(M)*noise_p+C_int_tmp);
            H_e=H_delay_obs(:,n0,t0);
            H_e=filter*H_e;
            He_delay_w(:,n0,t0)=H_e;
        end
    end
    H_delay_w=permute(He_delay_w,[3,1,2]); %% n_delay,M,NSlots
    NMSE_vector=zeros(NSlots,1);
    for i0=1:NSlots
        H_delay=H_delay_w(:,:,i0);
        H_est=Fm*H_delay;
        H_freq=Ht_freq_w(:,:,i0);
        err=H_freq-H_est;
        NMSE_vector(i0)=sum(abs(err(:)).^2)/sum(abs(H_freq(:)).^2);
    end
    nmse=mean(NMSE_vector(:));
end

