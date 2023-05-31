function freq_h = generate_freq_channel(C_cluster,cluster_freq_phase)
    M=size(C_cluster,1);
    [N,num_cluster]=size(cluster_freq_phase);
    freq_h=zeros(M,N);
    for c0=1:num_cluster
        C_tmp=C_cluster(:,:,c0);
        h_tmp=C_tmp^(1/2)*(sqrt(1/2).*(randn(M,1)+1j.*randn(M,1)));
        phase=cluster_freq_phase(:,c0);
        for n0=1:N
            freq_h(:,n0)=freq_h(:,n0)+h_tmp.*phase(n0);
        end
    end
    freq_h=freq_h.';
end
