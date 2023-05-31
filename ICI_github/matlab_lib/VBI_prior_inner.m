function [post] = VBI_prior_inner(obj_prior,post,channel,P0)
    %%
    aw=obj_prior.aw;
    bw=obj_prior.bw;
    bar_a=obj_prior.bar_a;
    bar_b=obj_prior.bar_b;
    pri=obj_prior.pri;
    qsi=post.qsi;
    %%
    aew=aw+qsi*2;
    q_h_square=sum(abs(channel.mu_h).^2,3)+channel.C_h.*2;
    bew=bw+q_h_square.*qsi./P0;
    bar_ae=bar_a+sum(1-qsi)*2;
    bar_be=bar_b+sum(q_h_square.*(1-qsi));
    qs_bar=psi(bar_ae)-log(bar_be)-bar_ae/bar_be*q_h_square+log(1-pri);
    qs_x=psi(aew)-log(bew)-log(P0)-aew./bew.*q_h_square./P0+log(pri)-qs_bar;
    qsi=1./(1+exp(-qs_x));
    %%
    post.aew=aew;
    post.bew=bew;
    post.bar_ae=bar_ae;
    post.bar_be=bar_be;
    post.qsi=qsi;
end

