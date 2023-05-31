function C_freq = freq_covariance_dmrs(C_cluster,cluster_freq_phase,pilot_idx)
    M=size(C_cluster,1);
    n_freq=length(pilot_idx);
    [N,num_cluster]=size(cluster_freq_phase);
    C_freq=zeros(M*n_freq,M*n_freq);
    for c0=1:num_cluster
        C_tmp=C_cluster(:,:,c0);
        for i0=1:n_freq
            ni=pilot_idx(i0);
            phasei=cluster_freq_phase(ni,c0);
            idx_i=(1:M)+M*(i0-1);
            for j0=1:n_freq
                nj=pilot_idx(j0); 
                idx_j=(1:M)+M*(j0-1);
                phasej=cluster_freq_phase(nj,c0);
                C_freq(idx_i,idx_j)=C_freq(idx_i,idx_j)+phasei*phasej'*C_tmp;
            end
        end
    end
%     C_freq=C_freq./N;
end
