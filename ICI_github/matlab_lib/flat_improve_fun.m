function [H_est_w] = flat_improve_fun(Rx_w,lambda,A)
    [M,ite]=size(Rx_w);
    At=A;
    %%
    H_est_w=zeros(M,ite);
    for a0=1:10
        yh_w=zeros(M,ite);
        for i0=1:ite
            Rx=Rx_w(:,i0);
            [yh] = cvx_lasso(Rx,lambda,At,M);
            yh_w(:,i0)=yh;
            H_est=At*yh;
            H_est_w(:,i0)=H_est;
        end        
        %%
        problem.M = unitaryfactory(M, 1);
        problem.cost  = @(x) pro_cost(x,Rx_w,yh_w,ite);
        problem.egrad = @(x) pro_egrad(x,Rx_w,yh_w,ite,M);      % notice the 'e' in 'egrad' for Euclidean
        options.verbosity=0;
        options.stopfun = @mystopfun;
%                 warning('off', 'manopt:getHessian:approx') ;
        [At, ~, ~, ~] = conjugategradient(problem,At,options);
%                 [At, ~, ~, ~] =trustregions(problem,At,options);
    end
end

function grad = pro_egrad(x,Rx_w,yh_w,ite,M)
    grad=zeros(M,M);
    for i0=1:ite
        grad=grad-Rx_w(:,i0)*yh_w(:,i0)'-yh_w(:,i0)*Rx_w(:,i0)';
    end
end

function cost = pro_cost(x,Rx_w,yh_w,ite)
    cost=0;
    for i0=1:ite
        cost=cost-2*real(Rx_w(:,i0)'*x*yh_w(:,i0));
    end
end

function [yh] = cvx_lasso(Rx,lambda,At,M)
     H_grid=At\Rx;
     cvx_begin quiet
     cvx_solver mosek
     variable yh(M) complex
        minimize( sum_square_abs(H_grid-yh)+lambda*norm(yh,1) )
     cvx_end
end
