function C_freq = freq_covariance(C_cluster,cluster_freq_phase)
    M=size(C_cluster,1);
    [N,num_cluster]=size(cluster_freq_phase);
    C_freq=zeros(M,M);
    for c0=1:num_cluster
        C_tmp=C_cluster(:,:,c0);
        for n0=1:N
            phase=cluster_freq_phase(n0,c0);
            C_freq=C_freq+phase*phase'*C_tmp;
        end
    end
    C_freq=C_freq./N;
end
