function cell_intef_cdl = interference_cell_channel(Intf_cell,cdl,delay_intf)
    num_user=Intf_cell.num_user;
    cell_intef_cdl= cell(1,num_user);
    for userIdx = 1:num_user
        cell_intef_cdl{userIdx} = copy_cdl(cdl);
        delay=delay_intf(userIdx);
        AoD=Intf_cell.AoD(userIdx);
        ZoD=Intf_cell.ZoD(userIdx);
        cell_intef_cdl{userIdx}.AnglesAoD=cell_intef_cdl{userIdx}.AnglesAoD+AoD;
        cell_intef_cdl{userIdx}.AnglesZoD=cell_intef_cdl{userIdx}.AnglesZoD+ZoD;
        cell_intef_cdl{userIdx}.PathDelays=cell_intef_cdl{userIdx}.PathDelays+delay;
        swapTransmitAndReceive(cell_intef_cdl{userIdx});
%         cell_intef_cdl{userIdx}.ReceiveAntennaArray.Element='isotropic';
%         cell_intef_cdl{userIdx}.TransmitAntennaArray.Element='38.901';
    end
end