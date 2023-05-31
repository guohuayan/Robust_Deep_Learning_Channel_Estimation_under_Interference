function [estChannelGrid,offset_self] = get_channel_snapshot_ofdm_track(waveformInfo,channel,carrier,offset_ref,is_home,wave_length)
%     T = double(waveformInfo.Nfft)*3;
    SampleRate=waveformInfo.SampleRate;
    T = wave_length*SampleRate;
    txWaveform = sqrt(1/2).*(rand(T,1)+1j.*rand(T,1));
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

