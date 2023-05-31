function seq = generate_ZC_pilot(N,u,v)
%% u: [0,1,...,29]
%% v: [0,1]
    p = prevprime(N);
    bar_q=p*(u+1)/31;
    q=floor(bar_q+0.5)+v*(-1)^(floor(2*bar_q));
    ZC = zadoffChuSeq(q,p);
    L=mod(N,p);
    seq=zeros(N,1);
    seq(1:p)=ZC;
    seq((p+1):end)=ZC(1:L);
end