function [xhat,xvar] = CBG_Est_MRF( X_Pri, Var_Pri, lambda,theta,sigVar)

% lambda = SigEst_R.lambda;%Sparsity
% theta = SigEst_R.theta;%Mean of the Gaussian distribution
% sigVar = SigEst_R.sigVar;%Variance of the Gaussian distribution


alpha = -abs(X_Pri - theta).^2 ./ (sigVar + Var_Pri) + abs(X_Pri).^2 ./ Var_Pri;

beta = lambda./(1-lambda) .* (Var_Pri./(sigVar + Var_Pri)) .* exp(alpha);%0.5* is removed for complex estiamtion

pii = 1./(1 + 1./beta);
gamma = (X_Pri./Var_Pri + theta./sigVar) ./ (1./Var_Pri + 1./sigVar);
nu = 1./(1./Var_Pri + 1./sigVar);

xhat = pii.*gamma;
xvar = (pii.*(nu + abs(gamma).^2) - pii.^2 .* abs(gamma).^2);

end