function estChannelGrid = get_channel_aftersync(channel,carrier,pathGains,sampleTimes,pathFilters,offset)
    estChannelGrid = nrPerfectChannelEstimate(carrier,pathGains,pathFilters,offset,sampleTimes);
    estChannelGrid=squeeze(estChannelGrid(:,end,:));
    estChannelGrid=estChannelGrid.*sqrt(prod(channel.ReceiveAntennaArray.Size));
end

