function [mu_pos,var_pos,pnz,var_pos_mean] = vamp_BG_Gauss(BG_mu,BG_var,pri_mu,pri_var,lambda)
    post_var_scale = BG_var+pri_var + eps;
    gamma = (pri_mu.*BG_var+BG_mu.*pri_var)./post_var_scale;
    nu = pri_var.*BG_var./post_var_scale;
%     nu=pri_var.*BG_var./(BG_var+pri_var);
%     gamma=(pri_mu./pri_var+BG_mu./BG_var).*nu;
    %%
    pnz = (1-lambda)./(lambda).*pri_var./post_var_scale.*exp((abs(pri_mu-BG_mu).^2./post_var_scale-abs(pri_mu).^2./pri_var));
    pnz = 1./(1+pnz);
    %%
%     tmp=divide_two_Gauss_pdf(pri_mu,BG_mu,pri_var+BG_var+eps,0,pri_var);
%     pnz=1./(1+1./(lambda./(1-lambda).*tmp));
    pnz(isnan(pnz)) = 0.999;
%     pz=pz/tmp;
%     pnz=pnz/tmp;  %% non-zero probability
    %%
    mu_pos=pnz.*gamma;
    var_pos=pnz.*(nu+real(gamma).^2)-abs(mu_pos).^2;
    var_pos=min(max(real(var_pos),1e-11),1e11);
    var_pos_mean=mean(var_pos);
end

% function f = complex_Gauss_pdf(x,mu,var)
%     tmp=x-mu;
%     f=exp(-abs(tmp)^2/var)/var/pi;
% end

function f = divide_two_Gauss_pdf(x,mu1,var1,mu2,var2)
    f=exp(-abs(x-mu1).^2./var1+abs(x-mu2).^2./var2)./var1.*var2;
end



