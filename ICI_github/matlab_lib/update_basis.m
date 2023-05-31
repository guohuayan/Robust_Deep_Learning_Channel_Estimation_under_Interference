function [At] = update_basis(Rx,H_w,Fh,At)
    [M,~]=size(Rx);
    %%      
    problem.M = unitaryfactory(M, 1);
    problem.cost  = @(x) pro_cost(x,Rx,H_w,Fh);
    problem.egrad = @(x) pro_egrad(x,Rx,H_w,Fh);      % notice the 'e' in 'egrad' for Euclidean
    options.verbosity=0;
    options.stopfun = @mystopfun;
%                 warning('off', 'manopt:getHessian:approx') ;
    [At, ~, ~, ~] = conjugategradient(problem,At,options);
%                 [At, ~, ~, ~] =trustregions(problem,At,options);
end

function grad = pro_egrad(x,Rx,H_w,Fh)
    grad=-H_w'*Fh'*Rx;
end

function cost = pro_cost(x,Rx,H_w,Fh)
    cost=-trace(x*Rx'*Fh*H_w);
end

