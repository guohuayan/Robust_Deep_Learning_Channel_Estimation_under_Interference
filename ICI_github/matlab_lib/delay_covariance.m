function C_delay = delay_covariance(C_cluster,cluster_freq_phase,support_tap)
    M=size(C_cluster,1);
    [N,num_cluster]=size(cluster_freq_phase);
    Fm=dftmtx(N)/sqrt(N);
    cluster_time_phase=Fm'*cluster_freq_phase;
    C_delay=zeros(M,M,length(support_tap));
    for c0=1:num_cluster
        C_tmp=C_cluster(:,:,c0);
        for n0=1:length(support_tap)
            idx=support_tap(n0);
            phase=cluster_time_phase(idx,c0);
            C_delay(:,:,n0)=C_delay(:,:,n0)+phase*phase'*C_tmp;
        end
    end
end
