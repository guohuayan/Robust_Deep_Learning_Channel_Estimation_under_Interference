function Yo = cvx_lasso_ici(H_obs,Cu)
%     [~,n]=size(A);
    [M,n_delay]=size(H_obs);
    [V,D] = eig(Cu);
    Ci=V*sqrt(inv(D));
    sigma2=trace(Cu)/M;
    lambda=real(1/sqrt(sigma2));
%     n=m1*m2;
%     z=reshape(Z,[n,1]);
    cvx_begin quiet
        variable Y(M,n_delay) complex
%         minimize( sum(square_abs(Ci*(H_obs-Y) ),'all' )+lambda*sum(abs(Y),'all') )
        minimize( sum(sum(square_abs(Ci*(H_obs-Y) )) )+lambda*sum(sum(abs(Y))) )
    cvx_end
    Yo=Y;
end