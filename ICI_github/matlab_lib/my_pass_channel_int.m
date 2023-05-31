function [rxeWaveform] = my_pass_channel_int(channel,txeWaveform,nRxAnts,power_int)
    rxeWaveform = channel(txeWaveform);
    rxeWaveform = rxeWaveform.*sqrt(nRxAnts).*sqrt(power_int);
end

