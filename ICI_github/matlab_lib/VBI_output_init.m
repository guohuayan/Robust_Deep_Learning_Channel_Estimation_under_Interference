function channel = VBI_output_init(CP_len,Mh)
    channel.mu_h_old=zeros(CP_len,Mh,2);
    channel.C_h_old=zeros(CP_len,Mh);
    channel.mu_h=channel.mu_h_old;
    channel.mu_ht=reshape(channel.mu_h,[CP_len,Mh*2]);
    channel.C_h=channel.C_h_old;
end

