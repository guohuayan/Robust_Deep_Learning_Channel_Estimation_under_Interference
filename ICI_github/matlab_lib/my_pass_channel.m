function [rxeWaveform,pathGains,sampleTimes] = my_pass_channel(channel,txeWaveform,nRxAnts)
    [rxeWaveform,pathGains,sampleTimes] = channel(txeWaveform);
    rxeWaveform = rxeWaveform.*sqrt(nRxAnts);
end

