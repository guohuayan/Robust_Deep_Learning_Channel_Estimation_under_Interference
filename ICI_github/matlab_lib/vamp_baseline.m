function xest = vamp_baseline(Y,Cw,lambda,BG_var,niter)
% support: the length of nonzero element

[n_delay, M] = size(Y);

% Cw=var*eye(M);

x_pir = zeros(n_delay, M);

BG_mu=zeros(n_delay*M,1);
x_var_pir=BG_var*lambda;

% rho=0.95;
% x_BG_pos_damp=zeros(N,1);
% x_var_ext_damp=0;
x_old=Y;
for t=1:niter
    
    [x_pos,x_var_pos] = LMMSE_module(Y,Cw,x_pir,x_var_pir);
    
    x_pos_vec=reshape(x_pos,[n_delay*M,1]);
    x_pir_vec=reshape(x_pir,[n_delay*M,1]);
    
    [x_ext,x_var_ext] = Gauss_pos_to_ext(x_pos_vec,x_var_pos,x_pir_vec,x_var_pir);
    
%     x_var_ext_damp=(1-rho)*x_var_ext_damp+rho*x_var_ext;
    
    [A_n,gamma_n,nu_n] = BG_para(x_ext,x_var_ext, lambda, BG_mu,BG_var);
    
    x_BG_pos=A_n.*gamma_n;
    x_BG_var=A_n.*(nu_n+abs(gamma_n).^2)-abs(A_n.*gamma_n).^2;
    
    x_BG_var=mean(x_BG_var);
    
    lambda=mean(A_n);
%     BG_var=mean(x_BG_var)/lambda;
%     x_BG_var=min(max(real(x_BG_var),1e-11),1e11);
    
%     x_BG_pos_damp=(1-rho).*x_BG_pos_damp+rho.*x_BG_pos;
    
%     [x_BG_pos,x_BG_var_tmp,x_BG_nzp,x_BG_var] = vamp_BG_Gauss(BG_mu,BG_var,x_ext,x_var_ext,lambda);
    
    [x_pir_vec,x_var_pir] = Gauss_pos_to_ext(x_BG_pos,x_BG_var,x_ext,x_var_ext);
    x_pir=reshape(x_pir_vec,[n_delay,M]);
    
	%%
    xest=reshape(x_BG_pos,[n_delay,M]);
    err=x_old-xest;
    err=sum(abs(err(:)).^2);
    if err<1e-3
        break
    end
    x_old=xest;
end

% xest=reshape(x_BG_pos,[n_delay,M]);
end


function [A_n,gamma_n,nu_n] = BG_para(Rhat, Rvar, lambda, theta, phi)

%Calcualte Problem dimensions
[N, T] = size(Rhat);

%Calculate posterior means and variances
post_var_scale = phi + Rvar + eps;
gamma_n = (Rhat.*phi+Rvar.*theta)./post_var_scale;
nu_n = Rvar.*phi./post_var_scale;

%Find posterior that the component x(n,t) is active
A_n = (1-lambda)./(lambda).*phi./nu_n.*exp((abs(Rhat-theta).^2./post_var_scale-abs(Rhat).^2./Rvar));
A_n = 1./(1+A_n);
A_n(isnan(A_n)) = 0.999;
% A_n(isnan(A_n)) = 1;
end

function [mu_ext,var_ext] = Gauss_pos_to_ext(mu_pos,var_pos,mu_pir,var_pir)
    a=mu_pos;A=var_pos;
    b=mu_pir;B=var_pir;
    var_ext=1/(1/A-1/B);
    var_ext=min(max(real(var_ext),1e-11),1e11);
    mu_ext=var_ext*(a./A-b./B);
end

function [mu_pos,var_pos] = LMMSE_module(Y,Cw,x_pir,x_var_pir)
    %% This module support the case when the noise var is not identical
    %% input:
    %% H is partial unitary matrix
    %% Cw is the noise covariance matrix
    %% mu_pir, var_pir: the priori mean and variance
    [n_delay, M]=size(Y);
    %%
    Cov=inv(1/x_var_pir*eye(M)+inv(Cw));
    mu_pos=zeros(n_delay, M);
    Y_obs=Y-x_pir;
    for n0=1:n_delay
        H_obs=reshape(Y_obs(n0,:),[M,1]);
%         Ch=x_var_pir*eye(M);       
        mu_pos(n0,:)=reshape(Cov/Cw*H_obs,[1,M]);
    end
    mu_pos=mu_pos+x_pir;
    var_pos=trace(Cov)/M;
    %%
%     Cw=diag(diag(Cw));  %% make sure it is diagonal matrix
%     tmp=(Cw+phi*x_var_pir*phi');
%     mu_pos=x_pir+x_var_pir*phi'/tmp*(y-phi*x_pir);
%     var_pos=inv(eye(N)/x_var_pir+phi'/Cw*phi);
%     var_pos=trace(var_pos)/N;
    var_pos=min(max(real(var_pos),1e-11),1e11);
end