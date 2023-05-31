function [carrier,pusch,harqProcesses,encodeULSCH,decodeULSCH] = system_setup(num_pilot)
    simParameters = [];             % Clear simParameters variable
    % simParameters.NFrames = 2;      % Number of 10ms frames
    % simParameters.SNRIn = [-5 0 5]; % SNR range (dB)
%     displaySimulationInformation = true;
%     perfectChannelEstimator = true;
    %% Carrier and PUSCH Configuration
    % Bandwidth, numerology (SCS and CP type) and other general parameters
    simParameters.NRB = 16;                % Bandwidth in number of resource blocks (52RBs at 15kHz SCS for 10MHz BW)
    simParameters.SubcarrierSpacing = 30;  % 15, 30, 60, 120, 240 (kHz)
    simParameters.CyclicPrefix = 'Normal'; % 'Normal' or 'Extended'
    simParameters.NCellID = 0;             % Cell identity
    simParameters.NTxAnts = 1;             % Number of transmit antennas
    simParameters.NRxAnts = 64;             % Number of receive antennas
    % UL-SCH/PUSCH parameters
    simParameters.PUSCH.TargetCodeRate = 308 / 1024;      % Code rate used to calculate transport block sizes
    simParameters.PUSCH.PRBSet = (0:simParameters.NRB-1); % PUSCH PRB allocation
    simParameters.PUSCH.SymbolSet = 0:13;            % PUSCH symbol allocation in each slot
    simParameters.PUSCH.NohPRB = 0;                  % Additional RE overhead per PRB
    simParameters.PUSCH.EnableHARQ = true;           % Enable/disable HARQ, if disabled, single transmission with RV=0, i.e. no retransmissions
    simParameters.PUSCH.Modulation = 'QPSK';         % 'pi/2-BPSK', 'QPSK', '16QAM', '64QAM', '256QAM'
    simParameters.PUSCH.NLayers = 1;                 % Number of PUSCH layers
    simParameters.PUSCH.RNTI = 1;                    % Radio Network Temporary Identifier
    simParameters.PUSCH.TransformPrecoding = false;  % Enable/disable transform precoding
    simParameters.PUSCH.TxScheme = 'nonCodebook';    % Transmission scheme ('nonCodebook','codebook')
    simParameters.PUSCH.NAntennaPorts = 1;           % Number of antenna ports for codebook based precoding
    simParameters.PUSCH.TPMI = 0;                    % Precoding matrix indicator for codebook based precoding
    % PUSCH DM-RS configuration
    simParameters.PUSCH.PUSCHMappingType = 'A';      % PUSCH mapping type ('A'(slot-wise),'B'(non slot-wise))
    simParameters.PUSCH.DMRSTypeAPosition = 2;       % Mapping type A only. First DM-RS symbol position (2,3)
    simParameters.PUSCH.DMRSLength = 1;              % Number of front-loaded DM-RS symbols (1(single symbol),2(double symbol))
    simParameters.PUSCH.DMRSAdditionalPosition = num_pilot-1;  % Additional DM-RS symbol positions (max range 0...3)
    simParameters.PUSCH.DMRSConfigurationType = 1;   % DM-RS configuration type (1,2)
    simParameters.PUSCH.NumCDMGroupsWithoutData = 2; % CDM groups without data
    simParameters.PUSCH.NIDNSCID = 0;                % Scrambling identity (0...65535)
    simParameters.PUSCH.NSCID = 0;                   % Scrambling initialization (0,1)
    simParameters.PUSCH.NRSID = 0;                   % Scrambling ID for low-PAPR sequences (0...1007)
    simParameters.PUSCH.GroupHopping = 'Disable';    % Hopping type ('Enable','Disable')
    % Define the propagation channel type
    simParameters.ChannelType = 'CDL'; % 'CDL' or 'TDL'
    %%
    carrier = nrCarrierConfig;
    carrier.SubcarrierSpacing = simParameters.SubcarrierSpacing;
    carrier.CyclicPrefix = simParameters.CyclicPrefix;
    carrier.NSizeGrid = simParameters.NRB;
    carrier.NCellID = simParameters.NCellID;
    pusch = simParameters.PUSCH;
    %%
    % The sample rate for the channel model is set using the value returned 
    % from <matlab:edit('nrOFDMInfo') nrOFDMInfo>.
%     waveformInfo = nrOFDMInfo(carrier);
%     channel.SampleRate = waveformInfo.SampleRate;
    %%
    if pusch.EnableHARQ
        % From PUSCH demodulation requirements in RAN WG4 meeting #88bis
        % (R4-1814062)
        rvSeq = [0 2 3 1];
    else
        % HARQ disabled - single transmission with RV=0, no retransmissions
        rvSeq = 0;
    end
    % Create UL-SCH encoder System object
    encodeULSCH = nrULSCH;
    encodeULSCH.MultipleHARQProcesses = true;
    encodeULSCH.TargetCodeRate = pusch.TargetCodeRate;
    % Create UL-SCH decoder System object
    % Use layered belief propagation for LDPC decoding, with half the number of
    % iterations as compared to the default for belief propagation decoding
    decodeULSCH = nrULSCHDecoder;
    decodeULSCH.MultipleHARQProcesses = true;
    decodeULSCH.TargetCodeRate = pusch.TargetCodeRate;
    decodeULSCH.LDPCDecodingAlgorithm = 'Layered belief propagation';
    decodeULSCH.MaximumLDPCIterationCount = 6; 
    %%
     % Specify the order in which we cycle through the HARQ processes
    NHARQProcesses = 16;
%     harqSequence = 1:NHARQProcesses;
%     % Initialize the state of all HARQ processes and reset the UL-SCH 
%     % decoder
    harqProcesses = hNewHARQProcesses(NHARQProcesses,rvSeq,1);
end

