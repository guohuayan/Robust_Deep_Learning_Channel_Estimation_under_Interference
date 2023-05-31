function [pathGains,sampleTimes,pathFilters,offset_self] = get_channel_sync(waveformInfo,channel)
    T = double(waveformInfo.Nfft)*12;
    txWaveform = sqrt(1/2).*(rand(T,1)+1j.*rand(T,1));
    [~,pathGains,sampleTimes] = channel(txWaveform);
    pathFilters = getPathFilters(channel);
    [offset_self,~] = nrPerfectTimingEstimate(pathGains,pathFilters);
end