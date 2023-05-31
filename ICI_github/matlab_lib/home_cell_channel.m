function home_cdl = home_cell_channel(Home_cell,cdl)
    home_cdl=copy_cdl(cdl);
    home_cdl.AnglesAoD=home_cdl.AnglesAoD+Home_cell.AoD;
    home_cdl.AnglesZoD=home_cdl.AnglesZoD+Home_cell.ZoD;
    swapTransmitAndReceive(home_cdl);
%     home_cdl.ReceiveAntennaArray.Element='isotropic';
%     home_cdl.TransmitAntennaArray.Element='38.901';
end