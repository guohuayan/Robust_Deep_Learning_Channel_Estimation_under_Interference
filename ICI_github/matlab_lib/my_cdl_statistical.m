function [cluster_freq_phase,C_cluster] = my_cdl_statistical(home_cdl,N,CP_len,CP_shift)
    home_cdlInfo = info(home_cdl);
    Fm=dftmtx(N)/sqrt(N);
    Mr=home_cdl.ReceiveAntennaArray.Size(2);
    Nr=home_cdl.ReceiveAntennaArray.Size(1);
    M=Mr*Nr;
    cASA=home_cdl.AngleSpreads(2);
    cZSA=home_cdl.AngleSpreads(4);
    gain=home_cdl.AveragePathGains;
    Pn=10.^(gain/10);
    num_cluster=length(Pn);
    Pn=Pn./sum(Pn);
    starting=home_cdlInfo.ChannelFilterDelay+round(CP_len*CP_shift);
    delay=home_cdl.PathDelays.*home_cdl.SampleRate+starting+1;
    cluster_freq_phase=zeros(N,num_cluster);
    for c0=1:num_cluster
        alpha=delay(c0)-1;
        phase=exp(-1j.*2*pi*(0:(N-1)).'/N*alpha)./sqrt(N);
        phase=Fm'*phase;
        [phase] = eliminate_delay_cluster(phase);
        phase=Fm*phase;
        cluster_freq_phase(:,c0)=phase;
    end
    cluster_freq_phase=cluster_freq_phase.*sqrt(N);
    alpha_w=[0.0447,0.1413,0.2492,0.3715,0.5129,0.6797,0.8844,1.1481,1.5195,2.1551];
    alpha_w=[alpha_w,-alpha_w];
    C_cluster=zeros(M,M,num_cluster);
    for c0=1:num_cluster
        alpha_wr=alpha_w;
        if (strcmpi(home_cdlInfo.ClusterTypes{c0},'LOS'))
            alpha_wr(:) = 0; % no ray offseting for LOS path (only first ray will be used anyway)
        end
        phi_AoA=alpha_wr.*cASA+home_cdl.AnglesAoA(c0);
        theta_ZoA=alpha_wr.*cZSA+home_cdl.AnglesZoA(c0);
        phi_AoA = wrapAzimuthAngles(phi_AoA); % wrap azimuth angles to [-180,180]
        theta_ZoA = wrapZenithAngles(theta_ZoA);
        Cg_c=zeros(M,M);
        power_ray=Pn(c0);
        power_ray=power_ray/length(alpha_wr);
        for r0=1:length(alpha_wr)    
            phi_AoD_r=phi_AoA(r0);
            theta_ZoD_r=theta_ZoA(r0);
            alpha_R=upa_vec(Mr,Nr,phi_AoD_r,theta_ZoD_r);
            C_tmp=power_ray.*alpha_R*alpha_R';
            Cg_c=Cg_c+C_tmp;
        end
        Cg_c=regular_smdf(Cg_c);
        C_cluster(:,:,c0)=Cg_c;
    end
end