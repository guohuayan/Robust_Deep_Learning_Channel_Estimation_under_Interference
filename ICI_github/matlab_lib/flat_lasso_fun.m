function [yh] = flat_lasso_fun(Rx,lambda,At,M)
     H_grid=At\Rx;
     cvx_begin quiet
     cvx_solver mosek
     variable yh(M) complex
        minimize( sum_square_abs(H_grid-yh)+lambda*norm(yh,1) )
     cvx_end
end

