function home_cdl = channel_3gpp(Info,cdl_in)
    home_cdl = clone(cdl_in);
%     cdlInfo = info(cdl_in);
    home_cdl.DelayProfile='Custom';
    home_cdl.AngleSpreads=Info.AngleSpreads;
    home_cdl.XPR=Info.XPR;
%     home_cdl.HasLOSCluster=Info.HasLOSCluster;
%     home_cdl.KFactorFirstCluster=Info.KFactorFirstCluster;
    home_cdl.PathDelays=Info.PathDelays;
    home_cdl.AveragePathGains=Info.AveragePathGains;
    home_cdl.AnglesAoD=Info.AnglesAoD;
    home_cdl.AnglesAoA=Info.AnglesAoA;
    home_cdl.AnglesZoD=Info.AnglesZoD;
    home_cdl.AnglesZoA=Info.AnglesZoA;
    %%
    swapTransmitAndReceive(home_cdl);
%     home_cdl.ReceiveAntennaArray.Element='isotropic';
%     home_cdl.TransmitAntennaArray.Element='38.901';
end