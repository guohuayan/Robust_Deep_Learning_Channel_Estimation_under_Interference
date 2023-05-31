function [rxeWaveform,pathGains,sampleTimes] = my_pass_channel_int_v2(channel,txeWaveform,nRxAnts,power_int)
    [rxeWaveform,pathGains,sampleTimes] = channel(txeWaveform);
    rxeWaveform = rxeWaveform.*sqrt(nRxAnts).*sqrt(power_int);
end

