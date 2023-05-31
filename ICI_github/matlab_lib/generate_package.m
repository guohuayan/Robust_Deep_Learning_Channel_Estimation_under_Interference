function [puschIndices,dmrsIndices,puschIndicesInfo,puschSymbols,puschGrid] = generate_package(carrier,pusch,harqProcesses,encodeULSCH)
    % Create 'ue' structure
    ue = struct();
    ue.NRB = carrier.NSizeGrid;
    ue.CyclicPrefix = carrier.CyclicPrefix;
    ue.SubcarrierSpacing = carrier.SubcarrierSpacing;
    ue.NCellID = carrier.NCellID;    
    % Timing offset, updated in every slot for perfect synchronization and
    % when the correlation is strong for practical synchronization
    offset = 0;
    nslot=1;
    % Update the carrier and PUSCH slot numbers to account for a new
    % PUSCH transmission
    carrier.NSlot = nslot;
    pusch.NSlot = nslot;

    % Calculate the transport block size for this slot
    [puschIndices,dmrsIndices,dmrsSymbols,puschIndicesInfo] = hPUSCHResources(ue,pusch);
    MRB = numel(pusch.PRBSet);
    TBS = nrTBS(pusch.Modulation,pusch.NLayers,MRB,puschIndicesInfo.NREPerPRB,pusch.TargetCodeRate,pusch.NohPRB);

    % Get HARQ process index for the current PUSCH from the HARQ index
    % table
    harqProcIdx=1;

    trBlk = randi([0 1],TBS,1);
    setTransportBlock(encodeULSCH,trBlk,harqProcIdx-1);

    % UL-SCH encoding
    codedTrBlock = encodeULSCH(pusch.Modulation,pusch.NLayers,puschIndicesInfo.G,harqProcesses(harqProcIdx).RV,harqProcIdx-1);

    % PUSCH modulation, including codebook based MIMO precoding if 
    % TxScheme = 'codebook'
    puschSymbols = nrPUSCH(codedTrBlock,pusch.Modulation,pusch.NLayers,carrier.NCellID,pusch.RNTI,pusch.TransformPrecoding,MRB,pusch.TxScheme,pusch.NAntennaPorts,pusch.TPMI);

    % Create resource grid associated with PUSCH transmission period
    nTxAnts=1;
    puschGrid = nrResourceGrid(carrier,nTxAnts);
end

