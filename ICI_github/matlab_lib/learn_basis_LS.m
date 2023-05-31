function A = learn_basis_LS(H_est)
    [~,M]=size(H_est);
    Mt=M/2;
    H1=H_est(:,1:Mt).';
    H2=H_est(:,(Mt+1):end).';
    Ct=H1*H1'+H2*H2';
    [A,~] = eig(Ct);
end

