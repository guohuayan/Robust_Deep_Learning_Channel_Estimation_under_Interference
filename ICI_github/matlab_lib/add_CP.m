function signal = add_CP(data,CP_len)
    [N,Nt]=size(data);
    signal=zeros(N+CP_len,Nt);
    signal(1:CP_len,:)=data((end-CP_len+1):end,:);
    signal((CP_len+1):end,:)=data;
end