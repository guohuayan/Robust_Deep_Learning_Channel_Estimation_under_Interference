function [At] = flat_update_basis(Dt,At)
    [M,~]=size(At);
    %%      
    problem.M = unitaryfactory(M, 1);
    problem.cost  = @(x) pro_cost(x,Dt);
    problem.egrad = @(x) pro_egrad(x,Dt);      % notice the 'e' in 'egrad' for Euclidean
    options.verbosity=0;
    options.stopfun = @mystopfun;
%                 warning('off', 'manopt:getHessian:approx') ;
    [At, ~, ~, ~] = conjugategradient(problem,At,options);
%                 [At, ~, ~, ~] =trustregions(problem,At,options);
end

function grad = pro_egrad(x,Dt)
    grad=-Dt';
end

function cost = pro_cost(x,Dt)
    cost=-trace(Dt*x);
end

