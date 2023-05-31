function data = remove_CP(signal,CP_len,shift)
    [L,~]=size(signal);
    N=L-CP_len;
%     data=zeros(N,Nr);
    shift_n=round(shift*CP_len);
    data=signal((CP_len-shift_n+1):(CP_len-shift_n+N),:);
end