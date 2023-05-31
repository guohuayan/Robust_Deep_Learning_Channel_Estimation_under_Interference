function cdl = copy_cdl(cdl_in)
    cdl = clone(cdl_in);
    cdlInfo = info(cdl_in);
    cdl.DelayProfile='Custom';
    if isequal(cdl_in.DelayProfile,'CDL-A')
        cdl.AngleSpreads =[5.0 11.0 3.0 3.0];
        cdl.XPR=10;
    elseif isequal(cdl_in.DelayProfile,'CDL-B')
        cdl.AngleSpreads =[10.0 22.0 3.0 7.0];
        cdl.XPR=8;
    elseif isequal(cdl_in.DelayProfile,'CDL-C')
        cdl.AngleSpreads =[2.0 15.0 3.0 7.0];
        cdl.XPR=7;
    elseif isequal(cdl_in.DelayProfile,'CDL-D')
        cdl.AngleSpreads =[5.0 8.0 3.0 3.0];
        cdl.XPR=11;
    elseif isequal(cdl_in.DelayProfile,'CDL-E')
        cdl.AngleSpreads =[5.0 11.0 3.0 7.0];
        cdl.XPR=8;
    end
    if isequal(cdlInfo.ClusterTypes{1},'LOS')
        cdl.HasLOSCluster=1;
        cdl.KFactorFirstCluster=cdlInfo.KFactorFirstCluster;
        cdlInfo.PathDelays(1)=[];
        cdlInfo.AveragePathGains(1)=[];
        cdlInfo.AnglesAoD(1)=[];
        cdlInfo.AnglesAoA(1)=[];
        cdlInfo.AnglesZoD(1)=[];
        cdlInfo.AnglesZoA(1)=[];
    end
    cdl.PathDelays=cdlInfo.PathDelays;
    cdl.AveragePathGains=cdlInfo.AveragePathGains;
    cdl.AnglesAoD=cdlInfo.AnglesAoD;
    cdl.AnglesAoA=cdlInfo.AnglesAoA;
    cdl.AnglesZoD=cdlInfo.AnglesZoD;
    cdl.AnglesZoA=cdlInfo.AnglesZoA;
%     cdl.TransmitAntennaArray.Element='isotropic';
%     cdl.ReceiveAntennaArray.Element='isotropic';
end