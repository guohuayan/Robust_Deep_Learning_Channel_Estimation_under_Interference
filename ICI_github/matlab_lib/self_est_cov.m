function Ch = self_est_cov(H_est)
    [N,M]=size(H_est);
    Mt=M/2;
    H1=H_est(:,1:Mt).';
    H2=H_est(:,(Mt+1):end).';
    Ct=H1*H1'+H2*H2';
    Ct=Ct./N./2;
    Ch=kron(eye(2),Ct);
end

