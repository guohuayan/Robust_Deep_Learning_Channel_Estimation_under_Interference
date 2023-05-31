function Ch = learn_basis_post(obj_prior,post,channel,P0)
    %%
    aw=obj_prior.aw;
    bw=obj_prior.bw;
    bar_a=obj_prior.bar_a;
    bar_b=obj_prior.bar_b;
    pri=obj_prior.pri;
    qsi=post.qsi;
    %%
    aew=qsi*2;
%     q_h_square=sum(abs(channel.mu_h).^2,3)+channel.C_h.*2;
    mu_h=channel.mu_ht;
%     Ch=channel.C_h;
%     [CP_len,Mh]=size(Ch);
%     qsi=post.qsi;
    h1=mu_h(:,1:Mh);
    h2=mu_h(:,(Mh+1):end);
%     %%
%     Ch=Ch.*qsi;
%     h1=h1.*qsi;
%     h2=h2.*qsi;
%     %%
    Ch=channel.C_h.*2;
%     Ch=real(diag(Ch));
    tmp=zeros(Mh,Mh);
    for i0=1:CP_len
        h1t=h1(i0,:).';
        h2t=h2(i0,:).';
        tmp=tmp+(h1t*h1t'+h2t*h2t');
    end
%     Ch=Ch+tmp/sum(mean(qsi,2));
%     Ch=kron(eye(2),Ch);
end