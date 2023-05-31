function [C_delay,power_tap] = delay_covariance_dmrs(C_cluster,cluster_freq_phase,Fm)
    M=size(C_cluster,1);
    num_cluster=size(C_cluster,3);
    [~,Len_delay]=size(Fm);
    cluster_time_phase=Fm'*cluster_freq_phase;
    C_delay=zeros(M,M,Len_delay);
    for c0=1:num_cluster
        C_tmp=C_cluster(:,:,c0);
        for n0=1:Len_delay
            phase=cluster_time_phase(n0,c0);
            C_delay(:,:,n0)=C_delay(:,:,n0)+phase*phase'*C_tmp;
        end
    end
    power_tap=zeros(Len_delay,1);
    for n0=1:Len_delay
        power_tap(n0)=real(trace(C_delay(:,:,n0)));
    end
end
