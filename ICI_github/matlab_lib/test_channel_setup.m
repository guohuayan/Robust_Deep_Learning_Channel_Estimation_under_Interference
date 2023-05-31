function [channel,Intercell_channel_G,power_int,Doppler] = test_channel_setup(speed,carrier)
    fc=3.5e9;
    waveformInfo = nrOFDMInfo(carrier);
    load('sim_multi_cell_interference.mat','G_Intf_cell','Home_cell');
    cdl = nrCDLChannel;
    cdl.DelayProfile = 'CDL-B';
    cdl.TransmitAntennaArray.Size = [4 8 2 1 1];
    cdl.ReceiveAntennaArray.Size = [1 1 1 1 1];
    cdl.TransmitAntennaArray.Element='38.901';
    cdl.ReceiveAntennaArray.Element='isotropic';
    cdlInfo = info(cdl);
    %%
    cspeed=physconst('LightSpeed');
%     speed=3;
    move=speed*1e3/3600;
    Doppler=(move/cspeed)*fc;
%     Doppler/30e3
    %% home_cell
    channel = home_cell_channel(Home_cell,cdl);
    channel.MaximumDopplerShift = Doppler; % in Hz
    channel.Seed=75;
    channel.CarrierFrequency=fc;
    channel.SampleRate = waveformInfo.SampleRate;
    %%
    Ci=size(G_Intf_cell,1);
    Intercell_channel_G=cell(Ci,1);
    cv=physconst('LightSpeed');
    %%
    for c0=1:Ci
        Intf_cell=G_Intf_cell{c0};
        distance_diff=0;
        delay_intf=distance_diff./cv;
    %     fprintf('start cell %d\n',c0);
        cell_intef_cdl_w = interference_cell_channel(Intf_cell,cdl,delay_intf);
        %% channel
        channel_i=cell_intef_cdl_w{1};
        channel_i.MaximumDopplerShift = Doppler; % in Hz
        channel.Seed=69+c0*2;
        channel_i.CarrierFrequency=fc;
        channel_i.SampleRate = waveformInfo.SampleRate;
        Intercell_channel_G{c0}=channel_i;
    end
    %% fix support and normlized RX sum power
    power_int=1/Ci;
end

