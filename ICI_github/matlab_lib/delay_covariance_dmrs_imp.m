function C_delay = delay_covariance_dmrs_imp(C_cluster,cluster_freq_phase,Fm)
    M=size(C_cluster,1);
    num_cluster=size(C_cluster,3);
    [~,Len_delay]=size(Fm);
    cluster_time_phase=Fm'*cluster_freq_phase;
    C_delay=zeros(M*Len_delay,M*Len_delay);
    for c0=1:num_cluster
        C_tmp=C_cluster(:,:,c0);
        for i0=1:Len_delay
            idx_i=(1:M)+M*(i0-1);
            phasei=cluster_time_phase(i0,c0);
            for j0=1:Len_delay
                idx_j=(1:M)+M*(j0-1);
                if abs(j0-i0)<=1
                    phasej=cluster_time_phase(j0,c0);
                else
                    phasej=0;
                end
                C_delay(idx_i,idx_j)=C_delay(idx_i,idx_j)+phasei*phasej'*C_tmp;
            end
        end
    end
end
