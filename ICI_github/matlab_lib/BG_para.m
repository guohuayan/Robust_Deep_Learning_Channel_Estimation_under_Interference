function [A_n,gamma_n,nu_n] = BG_para(x1_ext, x1_var_ext, lambda1, BG_mu, BG_var)

%Calcualte Problem dimensions
N = length(x1_ext);
Nv= length(x1_var_ext);
if Nv==1
    x1_var_ext=x1_var_ext.*ones(N,1);
end
Nv= length(BG_var);
if Nv==1
    BG_var=BG_var.*ones(N,1);
end
Nv= length(lambda1);
if Nv==1
    lambda1=lambda1.*ones(N,1);
end
%%

%Calculate posterior means and variances
post_var_scale = BG_var + x1_var_ext + eps;
gamma_n = (x1_ext.*BG_var+x1_var_ext.*BG_mu)./post_var_scale;
nu_n = x1_var_ext.*BG_var./post_var_scale;

%Find posterior that the component x(n,t) is active
% t1=abs(x_ext-BG_mu).^2./post_var_scale;
% t2=abs(x_ext).^2./x_var_ext;
A_n = (1-lambda1)./(lambda1).*BG_var./nu_n.*exp(abs(x1_ext-BG_mu).^2./post_var_scale-abs(x1_ext).^2./x1_var_ext);
A_n = 1./(1+A_n);
if sum(isnan(A_n))>0
    fprintf('nan of A_n\n');
%     pause
end
% A_n(isnan(A_n)) = 0.999;
end