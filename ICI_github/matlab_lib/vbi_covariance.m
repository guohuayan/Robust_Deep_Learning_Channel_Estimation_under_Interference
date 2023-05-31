function Ce = vbi_covariance(channel)
    [CP_len,M,group]=size(channel.mu_h);
    Ce=zeros(M,M);
    for i0=1:CP_len
        for j0=1:group
            h=channel.mu_h(i0,:,j0);
            h=reshape(h,[M,1]);
            Ce=Ce+h*h';
        end
    end
end

