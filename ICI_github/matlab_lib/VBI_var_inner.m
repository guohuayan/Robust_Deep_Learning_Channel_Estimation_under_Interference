function [var] = VBI_var_inner(post,P0)
    %%
    var.alpha_bar=(post.bar_ae/post.bar_be);
    var.alpha_w=(post.aew./post.bew);
    var.An_w=(post.qsi.*var.alpha_w/P0+(1-post.qsi)*var.alpha_bar);
    var.An_w=real(var.An_w);
end

