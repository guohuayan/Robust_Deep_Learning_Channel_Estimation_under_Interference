function seq_out = shift_ZC_pilot(seq,Nt,i)
%% 0=<i<Nt
    N=length(seq);
    C=12;
    n_cs=mod(C/Nt*i,C);
    alpha=2*pi*n_cs/C;
    shift=exp(-1j.*alpha.*(1:N));
    shift=reshape(shift,size(seq));
    seq_out=shift.*seq;
end