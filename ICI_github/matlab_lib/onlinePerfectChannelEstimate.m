function TrueChannel = onlinePerfectChannelEstimate(carrier,pathGains,pathFilters,offset,sampleTimes,nRxAnts)
    TrueChannel =nrPerfectChannelEstimate(carrier,pathGains,pathFilters,offset,sampleTimes).*sqrt(nRxAnts);
end

