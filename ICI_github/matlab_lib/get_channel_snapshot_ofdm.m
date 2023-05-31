function [estChannelGrid,offset_self] = get_channel_snapshot_ofdm(waveformInfo,channel,carrier,offset_ref,is_home)
    T = double(waveformInfo.Nfft)*3;
    txWaveform = sqrt(1/2).*(randn(T,1)+1j.*randn(T,1));
    [~,pathGains,sampleTimes] = channel(txWaveform);
    pathFilters = getPathFilters(channel);
    [offset_self,~] = nrPerfectTimingEstimate(pathGains,pathFilters);
    if is_home==1
        offset=offset_self;
    else
        offset=offset_ref;
    end
    estChannelGrid = myPerfectChannelEstimate(carrier,pathGains,pathFilters,offset,sampleTimes);
    estChannelGrid=squeeze(estChannelGrid(:,end,:));
    estChannelGrid=estChannelGrid.*sqrt(prod(channel.ReceiveAntennaArray.Size));
end

