function Yo = cvx_lasso_ici_free(Z,lambda)
%     [~,n]=size(A);
    [m1,m2]=size(Z);
    n=m1*m2;
    z=reshape(Z,[n,1]);
    cvx_begin quiet
        variable y(n) complex
%         minimize( sum_square_abs(z-A*y)+lambda*max(norm(y,1),p1) )
%         minimize( sum_square_abs(z-y)+lambda*norm(y,1) )
%         minimize( sum_square_abs(z-y)+lambda*sum(abs(y)) )
        minimize( sum(square_abs(z-y))+lambda*sum(abs(y)) )
    cvx_end
    yo=y;
    Yo=reshape(yo,[m1,m2]);
end