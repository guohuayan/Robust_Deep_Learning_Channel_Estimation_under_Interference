function flag = VBI_err_cal(channel,CP_len,M)
    err1=channel.mu_h_old-channel.mu_h;
    flag=(sum(abs(err1(:)).^2))/CP_len/M;
%     err2=channel.C_h_old-channel.C_h;
%     flag=(sum(abs(err1(:)).^2)+sum(abs(err2(:)).^2))/CP_len/M;
end

