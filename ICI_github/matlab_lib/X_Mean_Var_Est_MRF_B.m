function [theta_sig_update,var_sig_update,lambda_update] = X_Mean_Var_Est_MRF_B(S,X_Pri,Var_Pri,lambda,theta_sig,var_sig)
% lambda = SigEst_R.lambda;%Sparsity
% theta_sig = SigEst_R.theta; % Mean of the Gaussian distribution of x(channel coefficient)
% var_sig = SigEst_R.sigVar;%Variance of the Gaussian distribution of x

exp_factor_alpha = -abs(X_Pri - theta_sig).^2./((Var_Pri + var_sig)) + abs(X_Pri).^2./((Var_Pri));
exp_factor_beta = lambda./(1-lambda) .* (Var_Pri./(var_sig + Var_Pri)) .* exp(exp_factor_alpha);
% for ttt = 1:length(exp_factor_beta)
%     if(isnan(exp_factor_beta(ttt)))
%         exp_factor_beta(ttt) = 1e-9;
%     end
% end
Est_PI = 1./(1 + 1./exp_factor_beta);

Est_Gamma = (X_Pri./Var_Pri + theta_sig./var_sig)./(1./Var_Pri + 1./var_sig);

Est_Nu = 1./(1./Var_Pri + 1./var_sig);

theta_sig_update = sum(Est_PI.*Est_Gamma)/sum(Est_PI);

var_sig_update = sum(Est_PI.*(abs(theta_sig-Est_Gamma).^2+Est_Nu))/sum(Est_PI);

lambda_update = mean(Est_PI);
end