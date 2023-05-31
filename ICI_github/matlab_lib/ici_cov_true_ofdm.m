function C_int = ici_cov_true_ofdm(Rint)
    [Nd,M]=size(Rint);
%     len_He_int=size(He_int,1);
%     H_int_avg=mean(He_int.',2);
    Corr_int=zeros(M,M);
    for ci0=1:Nd
%         h_tmp=He_int(ci0,:).'-H_int_avg;
        h_tmp=Rint(ci0,:).';
        Corr_int=Corr_int+h_tmp*h_tmp';
    end
%     C_int=Corr_int./(len_He_int-1)+H_int_avg*H_int_avg';
%     C_int=Corr_int./(len_He_int)+H_int_avg*H_int_avg';
    C_int=Corr_int./Nd;
%     C_int=C_int*nfft/Np;        
%         C_int=C_int-eye(M)*noise_var;  %% There exists difference between
%         -noise_var at first, and later
%     [V_cint,D_cint]=eig(C_int);
%     D_cint=real(diag(D_cint));
%     D_cint=D_cint-noise_var;
%     idx=find(D_cint>0);
%     C_int=V_cint(:,idx)*diag(D_cint(idx))*V_cint(:,idx)';
%     C_int=C_int./trace(C_int).*p_org;
end

