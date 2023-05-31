function Info = get_para_3gpp(scenario,fc,d_2D,h_UT,h_BS)
    if isequal(scenario,'UMi')
        fc=max(fc,2);
        DS_mu_lg=-0.24*log10(1+ fc) - 6.83;
        DS_sig_lg=0.16*log10(1+ fc) + 0.28;
        ASD_mu_lg=-0.23*log10(1+ fc) + 1.53;
        ASD_sig_lg=0.11*log10(1+ fc) + 0.33;
        ASA_mu_lg=-0.08*log10(1+ fc) + 1.81;
        ASA_sig_lg=0.05*log10(1+ fc) + 0.3;
        ZSA_mu_lg=-0.04*log10(1+ fc) + 0.92;
        ZSA_sig_lg=-0.07*log10(1+ fc) + 0.41;
        XPR_mu=8.0;
        %%
        ZSD_mu_lg=max(-0.5, -3.1*(d_2D/1000)+ 0.01*max(h_UT-h_BS,0) +0.2);
        ZSD_sig_lg=0.35;
%         ZOD_offset=-10^(-1.5*log10(max(10, d_2D))+3.3);
        ZOD_offset=0;
        %%
        r_tau=2.1;
        xi=3;
        num_cluster=19;
        AngleSpreads =[10.0 22.0 3/8*10^ZSD_mu_lg 7.0];
    elseif isequal(scenario,'UMa')
        fc=max(fc,6);
        DS_mu_lg=-6.28 - 0.204* log10(fc);
        DS_sig_lg=0.39;
        ASD_mu_lg=1.5 - 0.1144* log10(fc);
        ASD_sig_lg=0.28;
        ASA_mu_lg=2.08 - 0.27* log10(fc);
        ASA_sig_lg=0.11;
        ZSA_mu_lg=-0.3236* log10(fc) + 1.512;
        ZSA_sig_lg=0.16;
        XPR_mu=7;
        %%
        ZSD_mu_lg=max(-0.5, -2.1*(d_2D/1000)-0.01*(h_UT - 1.5)+0.9);
        ZSD_sig_lg=0.49;
%         ZOD_offset=7.66*log10(fc)-5.96-10^( (0.208*log10(fc)- 0.782) * log10(max(25, d_2D))+(-0.13*log10(fc)+2.03) -0.07*(h_UT-1.5) );
        ZOD_offset=0;
        %%
        r_tau=2.3;
        xi=3;
        num_cluster=20;
        AngleSpreads =[2.0 15.0 3/8*10^ZSD_mu_lg 7.0];
    elseif isequal(scenario,'RMa')
        DS_mu_lg=-7.43;
        DS_sig_lg=0.48;
        ASD_mu_lg=0.95;
        ASD_sig_lg=0.45;
        ASA_mu_lg=1.52;
        ASA_sig_lg=0.13;
        ZSA_mu_lg=0.58;
        ZSA_sig_lg=0.37;
        XPR_mu=7;
        %%
        ZSD_mu_lg=max(-1, -0.19*(d_2D/1000)-0.01*(h_UT - 1.5) + 0.28);
        ZSD_sig_lg=0.30;
%         ZOD_offset=atand((35 - 3.5)/d_2D) -atand((35 - 1.5)/d_2D);
        ZOD_offset=0;
        %%
        r_tau=1.7;
        xi=3;
        num_cluster=10;
        AngleSpreads =[2.0 3.0 3/8*10^ZSD_mu_lg 3.0];
    elseif isequal(scenario,'InH')
        fc=max(fc,6);
        DS_mu_lg=-0.28* log10(1+fc) - 7.173;
        DS_sig_lg=0.10* log10(1+fc) + 0.055;
        ASD_mu_lg=1.62;
        ASD_sig_lg=0.25;
        ASA_mu_lg=-0.11* log10(1+fc) + 1.863;
        ASA_sig_lg=0.12* log10(1+fc) + 0.059;
        ZSA_mu_lg=-0.15* log10(1+fc) + 1.387;
        ZSA_sig_lg=-0.09* log10(1+fc) + 0.746;
        XPR_mu=10;
        %%
        ZSD_mu_lg=1.08;
        ZSD_sig_lg=0.36;
        ZOD_offset=0;
        %%
        r_tau=3;
        xi=3;
        num_cluster=19;
        AngleSpreads =[5.0 11.0 3/8*10^ZSD_mu_lg 9.0];
    end
    %%
    DS=once_realize_gauss_lg(DS_mu_lg,DS_sig_lg);
    ASD=once_realize_gauss_lg(ASD_mu_lg,ASD_sig_lg);
    ASA=once_realize_gauss_lg(ASA_mu_lg,ASA_sig_lg);
    ZSA=once_realize_gauss_lg(ZSA_mu_lg,ZSA_sig_lg);
    ZSD=once_realize_gauss_lg(ZSD_mu_lg,ZSD_sig_lg);
    ASD=min(ASD,104);
    ASA=min(ASA,104);
    ZSA=min(ZSA,52);
    ZSD=min(ZSD,52);
    XPR=XPR_mu;
    %%
    Info.HasLOSCluster=0;
    Info.KFactorFirstCluster=-inf;
    Info.ClusterTypes=cell(1,num_cluster);
    for n0=1:num_cluster
        Info.ClusterTypes{n0}='NLOS';
    end
    %%
    Info.XPR=XPR;
    Info.AngleSpreads=AngleSpreads;
    %%
    tau=-r_tau*DS*log(rand(1,num_cluster));
    tau=sort(tau,'ascend' );
    tau=tau-tau(1);
    Info.PathDelays=tau;
    %%
    Pn=exp(-tau.*(r_tau-1)./r_tau./DS).*10.^(-xi.*randn(1,num_cluster)./10);
    Pn=Pn./max(Pn);
    Info.AveragePathGains=10.*log10(Pn);
    %%
    AnglesAoA=2*ASA/1.4*sqrt(-log(Pn))./Cphi(num_cluster);
    Info.AnglesAoA=AnglesAoA.*random_BPSK(num_cluster)+bounded_randn(num_cluster)*ASA./7;
    AnglesAoD=2*ASD/1.4*sqrt(-log(Pn))./Cphi(num_cluster);
    Info.AnglesAoD=AnglesAoD.*random_BPSK(num_cluster)+bounded_randn(num_cluster)*ASD./7;
    %%
    AnglesZoA=-ZSA*log(Pn)./Ctheta(num_cluster);
    Info.AnglesZoA=AnglesZoA.*random_BPSK(num_cluster)+bounded_randn(num_cluster)*ZSA./7;
    AnglesZoD=-ZSD*log(Pn)./Ctheta(num_cluster);
    Info.AnglesZoD=AnglesZoD.*random_BPSK(num_cluster)+bounded_randn(num_cluster)*ZSD./7+ZOD_offset;
%     Info.AnglesZoD
    %%
    idx=Info.AveragePathGains>=-25;
    Info.AveragePathGains=Info.AveragePathGains(idx);
    Info.PathDelays=Info.PathDelays(idx);
    Info.AnglesAoA=Info.AnglesAoA(idx);
    Info.AnglesAoD=Info.AnglesAoD(idx);
    Info.AnglesZoA=Info.AnglesZoA(idx);
    Info.AnglesZoD=Info.AnglesZoD(idx);
    Info.ClusterTypes=Info.ClusterTypes(idx);
end

function [x] = Ctheta(num_cluster)
    num=[8 10 11 12 15 19 20 25];
    value=[0.889 0.957 1.031 1.104 1.1088 1.184 1.178 1.282];
    idx=num==num_cluster;
    x=value(idx);
end

function [x] = Cphi(num_cluster)
    num=[4 5 8 10 11 12 14 15 16 19 20 25];
    value=[0.779 0.860 1.018 1.090 1.123 1.146 1.190 1.211 1.226 1.273 1.289 1.358];
    idx=num==num_cluster;
    x=value(idx);
end

function [x] = bounded_randn(num_cluster)
    x=zeros(1,num_cluster);
    for n0=1:num_cluster
        while(1)
            tmp=randn(1,1);
            if abs(tmp)<3
                x(n0)=tmp;
                break
            end
        end
    end
end

function [x] = random_BPSK(len)
    x=randi([0,1],1,len);
    x=1-2.*x;
end

function [x] = once_realize_gauss_lg(mu,sigma)
    x_lg=sigma*bounded_randn(1)+mu;
    x=10^x_lg;
end

function [x] = once_realize_gauss(mu,sigma)
    x=sigma*bounded_randn(1)+mu;
end

