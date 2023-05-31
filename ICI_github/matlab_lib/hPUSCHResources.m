%hPUSCHResources 5G NR PUSCH, DM-RS, and PT-RS resource element indices, DM-RS and PT-RS values
%   [IND,DMRSIND,DMRS,PTRSIND,PTRS,INFO] = hPUSCHResources(UE,CHS) returns
%   the resource element (RE) indices for a 5G NR PUSCH transmission,
%   associated PUSCH DM-RS and PT-RS, given the time (OFDM symbols) and
%   frequency (PRBs) allocation of the PUSCH and the DM-RS configuration.
%   The 1-based linear PUSCH indices are returned in matrix IND. They are
%   defined relative to a three-dimensional RE grid representing a
%   14/12-symbol slot for NRB resource blocks (in the PUSCH numerology)
%   across the DM-RS ports/layers of the PUSCH. Each column of IND
%   represents the grid locations for a separate port or layer (the third
%   dimension of the RE grid). The DM-RS and PT-RS RE indices have the same
%   format and are returned in matrix DMRSIND and PTRSIND respectively. The
%   complex values of associated DM-RS sequence and PT-RS sequence are also
%   returned in matrix DMRS and PTRS respectively. Additional information
%   about DM-RS, PT-RS, and resourcing is returned in the structure INFO.
%
%   [IND,DMRSIND,DMRS,INFO] = hPUSCHResources(UE,CHS) returns the resource
%   element (RE) indices for a 5G PUSCH transmission, associated PUSCH
%   DM-RS symbols and indices. Additional information about DM-RS, PT-RS,
%   and resourcing is returned in the structure INFO.
%
%   The UE-wide settings input, UE, must be a structure including 
%   the fields:
%   NRB               - Number of resource blocks for the carrier
%                       or bandwidth part (in the PUSCH numerology)
%   CyclicPrefix      - Optional. Cyclic prefix length 
%                       ('Normal'(default),'Extended')
%   SubcarrierSpacing - Optional. Subcarrier Spacing (kHz)
%                       (15(default),30,60,120,240)
%   NCellID           - Optional. Physical layer cell identity (0...1007)
%                       (default 1)
%
%   The PUSCH specific input, CHS, must be a structure including the
%   fields:
%   NSlot                   - Optional. Slot number of PUSCH transmission
%                             (default 0, can be the absolute slot number)
%   PRBSet                  - PRBs allocated to the PUSCH (0-based indices)
%   PRBRefPoint             - Optional. The PRB reference point used for
%                             reference signal(s) generation in case of
%                             CP-OFDM. The nominal value is the start of
%                             bandwidth part in common resource grid
%                             (default 0)
%   SymbolSet               - OFDM symbols allocated to PUSCH within slot
%                             (including DM-RS symbols, 0-based indices,
%                             max range 0...13)
%   PortSet                 - DM-RS antenna ports used by PUSCH
%                             (0-based indices, max range 0...11, mapping
%                             to ports p=0...11 respectively)
%   NLayers                 - Optional. Number of transmission layers 
%                             (only used if PortSet was not defined and 
%                             then it selects ports=0...NLayers-1)
%   Modulation              - Modulation type of codeword
%                             ('pi/2-BPSK','BPSK','QPSK','16QAM','64QAM','256QAM')
%   TransformPrecoding      - Optional. Transform precoding for SC-FDMA
%                             (0(default),1)
%     Required when transform precoding is enabled:
%       NRSID               - Scrambling ID for low-PAPR sequences (0..1007)
%       GroupHopping        - Hopping type ('Enable','Disable')
%   TxScheme                - Optional. Transmission scheme
%                             ('NonCodebook'(default),'Codebook')
%     Required when transmission scheme is codebook based:
%       NAntennaPorts       - Number of antenna ports after codebook precoding
%       TPMI                - Transmitted precoding matrix indicator 
%                             (max range 0...27)
%   IntraSlotFreqHopping    - Optional. Intra-slot frequency hopping
%                             ('Disabled'(default),'Enabled')
%   InterSlotFreqHopping    - Optional. Inter-slot frequency hopping
%                             ('Disabled'(default),'Enabled')
%     Required when intra-slot/inter-slot frequency hopping is enabled:
%       RBOffset            - Second hop RB position index (0-based)
%   RNTI                    - Optional. Radio network temporary identifier
%                             (0...65535) (default 1)
%   DMRSConfigurationType   - DM-RS configuration type (1,2)
%                             (Transform precoded case always uses type 1)
%   NumCDMGroupsWithoutData - Optional. Number of CDM groups without data
%                             (1(default),2,3)
%                             (Transform precoded case always uses 2)
%   NIDNSCID                - DM-RS scrambling identity (0...65535)
%                             (only used when transform precoding is disabled)
%   NSCID                   - DM-RS scrambling initialization (0,1)
%                             (only used when transform precoding is disabled)
%   EnablePTRS              - Optional. Enable PT-RS (0(default),1)
%   PTRSTimeDensity         - Optional. Time density of PT-RS (L_PT-RS). It
%                             is one of {1,2,4} for CP-OFDM and is either 1
%                             or 2 for SC-FDMA (default 1)
%   PTRSFrequencyDensity    - Optional. Frequency density of PT-RS (K_PT-RS)
%                             (2(default),4) (Only used for CP-OFDM)
%   PTRSNumSamples          - Optional. Number of PT-RS samples in a group
%                             (NGroupSamp) (2(default),4)
%                             (Only used for SC-FDMA)
%   PTRSNumGroups           - Optional. Number of PT-RS groups (NPTRSGroup)
%                             (2(default),4,8) (Only used for SC-FDMA)
%                             (When PTRSNumGroups is set to 8,
%                             PTRSNumSamples must be 4)
%   PTRSREOffset            - Optional. Resource element offset
%                             ('00'(default),'01','10','11') (Only used for
%                             CP-OFDM)
%   PTRSPortSet             - Optional. Antenna ports for PT-RS. The value
%                             must be a subset of DM-RS antenna ports. It
%                             is a scalar or two-element vector with values
%                             in range 0 to 3 for DM-RS configuration type
%                             1 and in range 0 to 5 for DM-RS configuration
%                             type 2 (default is lowest DM-RS port number)
%                             (Only used for CP-OFDM)
%   PTRSNID                 - Optional. PT-RS scrambling identity (0...1007)
%                             (default is the DM-RS scrambling identity NRSID)
%                             (Only used for SC-FDMA)
%
%   The DM-RS OFDM symbol locations can either be defined explicitly 
%   using the single parameter:
%   DMRSSymbolSet           - OFDM symbols containing the DM-RS within the 
%                             PUSCH allocation in slot 
%                             (0-based indices, max range 0...13)
%   or, defined implicitly via the following group of DM-RS parameters:
%   PUSCHMappingType        - PUSCH mapping type
%                             ('A'(slot-wise),'B'(non slot-wise))
%   DMRSTypeAPosition       - Mapping type A only. First DM-RS symbol position
%                             (2,3)
%   DMRSLength              - Number of consecutive front-loaded DM-RS OFDM
%                             symbols (1(single symbol),2(double symbol))
%   DMRSAdditionalPosition  - Additional DM-RS symbol positions
%                             (max range 0...3)
%
%   In terms of frequency domain DM-RS density, there are two different RRC
%   signaled configuration types ('dmrs-Type'). Configuration type 1 
%   defines 6 subcarriers per PRB per antenna port, comprising alternating
%   subcarriers. Configuration type 2 defines 4 subcarriers per PRB per
%   antenna ports, consisting of 2 groups of 2 neighboring subcarriers.
%   Different shifts are applied to the sets of subcarriers used, dependent
%   on the associated antenna port or CDM group. For type 1, there are 
%   2 possible CDM groups/shifts across up to 8 possible antenna ports
%   (p=0...7), and, for type 2, there are 3 possible CDM groups/shifts
%   across 12 ports (p=0...11). See TS 38.211 section 6.4.1.1 for the full
%   configuration details.
%
%   In terms of the time-domain DM-RS symbol positions, the PUSCH mapping
%   type ('PUSCHMappingType') can be either slot-wise (type A) or non
%   slot-wise (type B). When a UE is scheduled to receive PUSCH by a DCI,
%   this mapping type is signaled by the time-domain resource field in the
%   grant. The field acts as an index into an RRC configured table where
%   each row in the table specifies a combination of mapping type, slot
%   offset, K0, the symbol start and length indicator, SLIV. The mapping
%   type specifies the relative locations of the associated DM-RS. For
%   slot-wise mapping type A, the first DM-RS symbol is signaled by a field
%   in the MIB to be either 2 or 3 ('dmrs-TypeA-Position'). For the non
%   slot-wise mapping type B, the first DM-RS symbol is always the first
%   symbol of the PUSCH time allocation.
% 
%   The maximum number of DM-RS OFDM symbols used by a UE is configured by
%   RRC signaling ('dmrs-AdditionalPosition' and 'maxLength'). The DM-RS
%   can be a set of single symbols, distributed roughly uniformly across
%   the allocated PUSCH symbols, or 1 or 2 pairs of neighboring or 'double
%   symbol' DM-RS. The 'maxLength' RRC parameter (1 or 2 respectively)
%   configures whether only single symbol DM-RS or either single or double
%   symbol DM-RS can be used. In the latter case, the actual selection is
%   signaled in the DCI format 0_1 message. The 'dmrs-AdditionalPosition'
%   higher-layer parameter defines the number of additional single or
%   double symbol DM-RS that are transmitted. The valid combinations of
%   these two parameters is given by TS 38.211 tables 6.4.1.1.3-3,
%   6.4.1.1.3-4 and 6.4.1.1.3-6. In this function, the value of the
%   'DMRSLength' input parameter directly controls whether either single or
%   double symbols are used.
% 
%   INFO is the output structure containing the fields:
%   G             - Bit capacity of the PUSCH. This should be the length
%                   of codeword to be output from the UL-SCH transport
%                   channel
%   Gd            - Number of resource elements per layer/port, equal to
%                   the number of rows in the PUSCH indices
%   NREPerPRB     - Number of RE per PRB allocated to PUSCH
%   DMRSSymbolSet - The symbol numbers in a slot containing DM-RS (0-based)
%   CDMGroups     - CDM groups associated with the DM-RS antenna ports
%   CDMLengths    - A 2-element row vector [FD TD] specifying the length of
%                   FD-CDM and TD-CDM despreading required during channel
%                   estimation. The values depend on the frequency and
%                   time masks applied to groups of antenna ports
%   PTRSSymbolSet - The symbol numbers in a slot containing PT-RS (0-based)
%
%   Example:
%   % Display the locations of the PUSCH and PUSCH DM-RS resource elements 
%   % in a slot. 
%     
%   % Set the number of uplink carrier or BWP resource blocks and the 
%   % numerology (subcarrier spacing and cyclic prefix)
%   ue = struct('NRB',50);
%   ue.SubcarrierSpacing = 15;
%   ue.CyclicPrefix = 'Normal';
%   
%   % Get the number of OFDM symbols in a slot
%   symbperslot = sum(strcmpi(ue.CyclicPrefix,["Normal","Extended"]) .* [14 12]);
%   
%   % Specify the basic PUSCH allocation properties to be full band, full
%   % slot using 2 layers/ports (no transform precoding or frequency hopping)
%   pusch = struct();
%   pusch.NSlot = 0;                   % Slot number
%   pusch.PRBSet = 0:ue.NRB-1;         % Full band PRBs allocated to the PUSCH
%   pusch.SymbolSet = 0:symbperslot-1; % Full slot allocation top the PUSCH
%   pusch.PortSet = [0,2];             % Use DM-RS ports p=0 and p=2
%   pusch.Modulation = 'QPSK';         % Modulation scheme
% 
%   % Configure the PUSCH DM-RS for config type 1 and slot-wise, type A
%   % mapping. Use double symbols for front-loaded DM-RS and configure an
%   % additional symbol pair towards the end of the slot
%   pusch.DMRSConfigurationType = 1;  % DM-RS configuration type (1,2)
%   pusch.NIDNSCID = 1;               % DM-RS scrambling identity (0...65535)
%   pusch.NSCID = 0;                  % DM-RS scrambling initialization (0,1)
%   pusch.PUSCHMappingType = 'A';     % Slot-wise PUSCH mapping type
%   pusch.DMRSTypeAPosition = 2;      % First DM-RS symbol position for type A
%   pusch.DMRSLength = 2;             % Specify double front-loaded DM-RS
%   pusch.DMRSAdditionalPosition = 1; % Specify an additional DM-RS pair
%   pusch.NumCDMGroupsWithoutData = 3;% CDM groups without data
%
%   % Configure the PT-RS with frequency density set to 2, time density set
%   % to 1, resource element offset set to '00' and PT-RS antenna port set
%   % to 0
%   pusch.EnablePTRS = 1;             % Enable or disable PT-RS (1 or 0)
%   pusch.PTRSFrequencyDensity = 2;   % Frequency density (2,4)
%   pusch.PTRSTimeDensity = 1;        % Time density (1,2,4)
%   pusch.PTRSREOffset = '00';        % Resource element offset ('00','01','10','11')
%   pusch.PTRSPortSet = 0;            % Antenna port set for PT-RS
%
%   % Display PUSCH, DM-RS and PT-RS RE locations of the first port of
%   % the grid
%   slotgrid = zeros(12*ue.NRB,symbperslot,length(pusch.PortSet));
%   [ind,dmrsind,dmrs,ptrsind,ptrs,info] = hPUSCHResources(ue,pusch);
%   slotgrid(ind) = 20;                 % Use light blue for PUSCH RE 
%   slotgrid(dmrsind) = 40*abs(dmrs);   % Use green for DM-RS RE
%   slotgrid(ptrsind) = 70*abs(ptrs);   % Use yellow for PT-RS RE
%   figure;
%   imagesc(slotgrid(:,:,1));
%   title('PUSCH, DM-RS and PT-RS resource elements');
%   axis('xy'), xlabel('Symbols'), ylabel('Subcarriers');

%   Copyright 2019-2020 The MathWorks, Inc.

%#codegen

function [puschIndices,dmrsIndices,dmrsSymbols,varargout] = hPUSCHResources(ue,pusch)
 
    % Argument check 
    narginchk(2,2);
 
    % Get the carrier and pusch objects from the input structures
    [carrier,puschObj,indNonContig] = convertStructToObjects(ue,pusch);
    numSymNonContig = numel(indNonContig); % Number of OFDM symbols that are not part of the contiguous symbol allocation

    % Get the number of antenna ports
    if strcmpi(puschObj.TransmissionScheme,'nonCodebook')
        nports = puschObj.NumLayers;
    else
        nports = puschObj.NumAntennaPorts;
    end

    % Get the PUSCH resource element (RE) indices for the contiguous symbol
    % allocation
    [puschIndicesAll,gInfo] = nrPUSCHIndices(carrier,puschObj);
    puschIndicesAllLayers = double(puschIndicesAll);

    % Get the PUSCH RE indices
    nREPerPRB = gInfo.NREPerPRB;
    if numSymNonContig && ~isempty(puschIndicesAllLayers)
        % Remove the RE indices of the OFDM symbol locations, which are not
        % part of allocated OFDM symbols
        opts = struct;
        opts.IndexBase = '1based';
        opts.IndexStyle = 'subscript';
        gridSize = [carrier.NSizeGrid*12 carrier.SymbolsPerSlot nports];
        indSub = nr5g.internal.applyIndicesOptions(gridSize,opts,puschIndicesAll(:,1)-1);
        lm = repmat(indSub(:,2)-1,1,numSymNonContig) == repmat(indNonContig',numel(puschIndicesAll(:,1)),1);
        puschIndicesTemp = puschIndicesAllLayers(~sum(lm,2),:);
        % Get the number of REs available for data in the starting resource
        % block (RB), excluding DM-RS and PT-RS (if any)
        indSubWithinSymbolSet = nr5g.internal.applyIndicesOptions(gridSize,opts,puschIndicesTemp(:,1)-1);
        nREPerPRB = nnz((indSubWithinSymbolSet(:,3) == 1) & (indSubWithinSymbolSet(:,1) <= 12));
    else
        puschIndicesTemp = puschIndicesAllLayers;
    end
    puschIndices = reshape(puschIndicesTemp,[],nports);

    % Get the DM-RS OFDM symbol locations within the OFDM symbols allocated
    % for PUSCH
    nDMRSSym = numel(gInfo.DMRSSymbolSet);
    if numSymNonContig
        % Remove the DM-RS OFDM symbol locations, which are not part of
        % allocated OFDM symbols
        dmrsLogicalMatrix = repmat(gInfo.DMRSSymbolSet,numSymNonContig,1) == repmat(indNonContig,1,nDMRSSym);
        dmrsLogicalIndex = ~sum(dmrsLogicalMatrix,1);
    else
        dmrsLogicalIndex = true(nDMRSSym,1);
    end
    dmrsSymbolSet = gInfo.DMRSSymbolSet(dmrsLogicalIndex);

    % Get the PUSCH DM-RS symbols and indices, for the OFDM symbols
    % carrying DM-RS within the allocated OFDM symbols
    if ~isempty(dmrsSymbolSet)
        dmrsIndicesTemp = reshape(double(nrPUSCHDMRSIndices(carrier,puschObj)),[],nDMRSSym,nports);
        dmrsSymbolsTemp = reshape(nrPUSCHDMRS(carrier,puschObj),[],nDMRSSym,nports);
        dmrsIndices = reshape(dmrsIndicesTemp(:,dmrsLogicalIndex,:),[],nports);
        dmrsSymbols = reshape(dmrsSymbolsTemp(:,dmrsLogicalIndex,:),[],nports);
    else
        dmrsIndices = zeros(0,nports);
        dmrsSymbols = zeros(0,nports);
    end

    % Capture the number of PT-RS antenna ports
    if puschObj.TransformPrecoding
        % When transform precoding is enabled, the number of PT-RS ports
        % is equal to the number of layers
        nPTRSPorts = puschObj.NumLayers;
    else
        % When transform precoding is disabled, the number of PT-RS ports
        % depends on both the transmission scheme and the PT-RS antenna
        % ports configured
        if strcmpi(puschObj.TransmissionScheme,'nonCodebook')
            % For non-codebook-based transmission scheme, the number of
            % PT-RS ports is equal to the number of PT-RS ports configured
            if ~isempty(puschObj.PTRS.PTRSPortSet)
                nPTRSPorts = numel(puschObj.PTRS.PTRSPortSet);
            else
                nPTRSPorts = 1;
            end
        else
            % For codebook-based transmission scheme, the number of PT-RS
            % ports is equal to the number of antenna ports
            nPTRSPorts = nports;
        end
    end

    % Get the PUSCH PT-RS symbols and indices
    nPTRS = 0;
    if ~isempty(gInfo.PTRSSymbolSet) && ~isempty(dmrsSymbolSet)
        % Get the PT-RS OFDM symbol locations within the allocated OFDM
        % symbols
        symbolSet = puschObj.SymbolAllocation(1):sum(puschObj.SymbolAllocation(:))-1;
        ptrsSymbolSetTemp = nr5g.internal.pxsch.ptrsSymbolIndicesCPOFDM(symbolSet(symbolSet < carrier.SymbolsPerSlot), dmrsSymbolSet, puschObj.PTRS.TimeDensity);
        nPTRSSym = numel(ptrsSymbolSetTemp);
        if numSymNonContig
            % Remove the PT-RS OFDM symbol locations, which are not part of
            % allocated OFDM symbols
            ptrsLogicalMatrix = repmat(ptrsSymbolSetTemp,numSymNonContig,1) == repmat(indNonContig,1,nPTRSSym);
            ptrsLogicalIndex = ~sum(ptrsLogicalMatrix,1);
            % Assign the PUSCH DM-RS OFDM symbol locations, which are
            % within the allocated OFDM symbols to CustomSymbolSet
            puschObj.DMRS.CustomSymbolSet = dmrsSymbolSet;
        else
            ptrsLogicalIndex = true(nPTRSSym,1);
        end
        ptrsSymbolSet = ptrsSymbolSetTemp(ptrsLogicalIndex);

        % Generate PT-RS symbols and indices for the contiguous symbol
        % allocation
        ptrsInd = double(nrPUSCHPTRSIndices(carrier,puschObj));
        ptrsSym = nrPUSCHPTRS(carrier,puschObj);
        ptrsIndicesTemp = reshape(ptrsInd,[],nPTRSSym,nPTRSPorts);
        ptrsSymbolsTemp = reshape(ptrsSym,[],nPTRSSym,nPTRSPorts);
        powerFactor = 1;
        if puschObj.TransformPrecoding
            % Get PT-RS power factor for DFT-s-OFDM, depending on
            % modulation scheme
            powerFactor = ptrsPowerFactorDFTsOFDM(puschObj.Modulation);
        end

        % Get the PUSCH PT-RS symbols and indices, for the OFDM symbols
        % carrying PT-RS within the allocated OFDM symbols
        ptrsIndices = reshape(ptrsIndicesTemp(:,ptrsLogicalIndex,:),[],nPTRSPorts);
        ptrsSymbols = reshape(powerFactor*ptrsSymbolsTemp(:,ptrsLogicalIndex,:),[],nPTRSPorts);
        if puschObj.TransformPrecoding
            nPTRS = size(ptrsIndices,1);
        end
    else
        ptrsSymbolSet = zeros(1,0);
        ptrsIndices = zeros(0,nPTRSPorts);
        ptrsSymbols = zeros(0,nPTRSPorts);
    end

    % When the OFDM symbols allocated for PUSCH are non-contiguous, get the
    % number of REs in a resource block, depending on values of both the
    % properties EnablePTRS and TransformPrecoding
    if numSymNonContig && puschObj.EnablePTRS && ~puschObj.TransformPrecoding
        % Resource block offset, kRBRef
        nPUSCHRB = numel(unique(puschObj.PRBSet(:)));
        if mod(nPUSCHRB,puschObj.PTRS.FrequencyDensity) == 0
            kRBRef = mod(puschObj.RNTI,puschObj.PTRS.FrequencyDensity);
        else
            kRBRef = mod(puschObj.RNTI,mod(nPUSCHRB,puschObj.PTRS.FrequencyDensity));
        end
        % Add the number of REs available for data with number of PT-RS
        % symbols, when PT-RS is present in starting RB
        % Only one resource element is used for PT-RS in an RB
        nREPerPRB = nREPerPRB + numel(ptrsSymbolSet)*(~kRBRef);
    end

    % Combine information into output structure
    Gd = size(puschIndices,1)-nPTRS;
    qm = nr5g.internal.getQm(puschObj.Modulation);
    puschInfo.G = Gd*puschObj.NumLayers*qm;
    puschInfo.Gd = Gd;
    puschInfo.NREPerPRB = nREPerPRB;
    puschInfo.DMRSSymbolSet = dmrsSymbolSet;
    puschInfo.CDMGroups = puschObj.DMRS.CDMGroups;
    puschInfo.CDMLengths = puschObj.DMRS.CDMLengths;
    puschInfo.PTRSSymbolSet = ptrsSymbolSet;

    % Assign the outputs based on number of output arguments
    if nargout <= 4
        % [puschIndices,dmrsIndices,dmrsSymbols,puschInfo]
        varargout{1} = puschInfo;
    else
        % [puschIndices,dmrsIndices,dmrsSymbols,ptrsIndices,ptrsSymbols,puschInfo]
        varargout{1} = ptrsIndices;
        varargout{2} = ptrsSymbols;
        varargout{3} = puschInfo;
    end
end

function beta = ptrsPowerFactorDFTsOFDM(modulation)
%ptrsPowerFactorDFTsOFDM Provides the power factor (betaPrime) of PT-RS
% for DFT-s-OFDM according to TS 38.214 Table 6.2.3.2-1, based on the input
% modulation

    switch lower(modulation)
        case '16qam'
            beta = 3/sqrt(5);
        case '64qam'
            beta = 7/sqrt(21);
        case '256qam'
            beta = 15/sqrt(85);
        otherwise
            beta = 1;
    end

end

function [carrier,pusch,ind] = convertStructToObjects(ue,chs)
%convertStructToObjects Provides the configuration objects for given structures
%
%   [CARRIER,PUSCH,IND] = convertStructToObjects(UE,CHS) provides the
%   carrier configuration object CARRIER and physical uplink shared channel
%   object PUSCH, given the input structures UE-wide settings UE and
%   channel specific transmission configuration CHS. This function also
%   provides the OFDM symbol indices IND, which are not part of the
%   contiguous symbol allocation.

    % Get the carrier configuration with the inputs ue and chs
    carrier = nrCarrierConfig;
    carrier.NSizeGrid = ue.NRB;
    carrier = passign(chs,carrier,'PRBRefPoint','NStartGrid');
    carrier = passign(ue,carrier,'NCellID');
    carrier = passign(ue,carrier,'SubcarrierSpacing');
    carrier = passign(ue,carrier,'CyclicPrefix');
    carrier = passign(chs,carrier,'NSlot');

    % Get the pusch configuration object with the chs input
    pusch = nrPUSCHConfig;
    dmrs = nrPUSCHDMRSConfig;
    dmrs = passign(chs,dmrs,'PortSet','DMRSPortSet');
    numDMRSPorts = numel(dmrs.DMRSPortSet);
    if numDMRSPorts
        % Assign the number of layers to the number of DM-RS antenna ports
        pusch.NumLayers = numDMRSPorts;
    else
        % Get the number of layers, when DM-RS antenna port set is empty
        pusch = passign(chs,pusch,'NLayers','NumLayers');
    end
    pusch.Modulation = chs.Modulation;
    pusch = passign(chs,pusch,'PUSCHMappingType','MappingType');
    % Get SymbolAllocation value, depending on the values of SymbolSet
    % field
    if ~isfield(chs,'SymbolSet')
        pusch.SymbolAllocation = [0 carrier.SymbolsPerSlot];
        ind = zeros(0,1);
    elseif isempty(chs.SymbolSet)
        pusch.SymbolAllocation = [0 0];
        ind = zeros(0,1);
    else
        [lb,ub] = bounds(chs.SymbolSet);
        pusch.SymbolAllocation = [lb ub-lb+1];
        symTemp = lb:ub;
        logicalMatrix = repmat(symTemp,numel(chs.SymbolSet),1) == repmat(chs.SymbolSet(:),1,pusch.SymbolAllocation(2));
        ind = symTemp(~sum(logicalMatrix,1))';
    end
    pusch.PRBSet = chs.PRBSet;
    pusch = passign(chs,pusch,'TransformPrecoding');
    pusch = passign(chs,pusch,'TxScheme','TransmissionScheme');
    pusch = passign(chs,pusch,'TPMI');
    pusch = passign(chs,pusch,'NAntennaPorts','NumAntennaPorts');
    % Get FrequencyHopping value, depending on the values of both
    % InterSlotFreqHopping and IntraSlotFreqHopping fields
    if isfield(chs,'InterSlotFreqHopping') && strcmpi(chs.InterSlotFreqHopping,'enabled')
        pusch.FrequencyHopping = 'interSlot';
    elseif isfield(chs,'IntraSlotFreqHopping') && strcmpi(chs.IntraSlotFreqHopping,'enabled')
        pusch.FrequencyHopping = 'intraSlot';
    end
    pusch = passign(chs,pusch,'RBOffset','SecondHopStartPRB');
    pusch = passign(chs,pusch,'RNTI');

    % Set DM-RS parameters
    dmrs = passign(chs,dmrs,'DMRSConfigurationType','DMRSConfigurationType',~pusch.TransformPrecoding);
    dmrs = passign(chs,dmrs,'DMRSTypeAPosition');
    dmrs = passign(chs,dmrs,'DMRSAdditionalPosition');
    dmrs = passign(chs,dmrs,'DMRSLength');
    dmrs = passign(chs,dmrs,'DMRSSymbolSet','CustomSymbolSet');
    [dmrs,fieldPresent] = passign(chs,dmrs,'NIDNSCID','NIDNSCID',~pusch.TransformPrecoding);
    if ~fieldPresent
        % Assign NIDNSCID with physical layer cell identity NCellID, when
        % NIDNSCID field is not present or is empty
        dmrs.NIDNSCID = carrier.NCellID;
    end
    dmrs = passign(chs,dmrs,'NSCID','NSCID',~pusch.TransformPrecoding);
    % Get GroupHopping and SequenceHopping property values, depending on
    % GroupHopping field
    if isfield(chs,'GroupHopping') && pusch.TransformPrecoding
        if strcmpi(chs.GroupHopping,'enable')
            dmrs.GroupHopping = 1;
        elseif strcmpi(chs.GroupHopping,'disable')
            dmrs.SequenceHopping = 1;
        end
    end
    [dmrs,fieldPresent] = passign(chs,dmrs,'NRSID','NRSID',pusch.TransformPrecoding);
    if ~fieldPresent
        % Assign NRSID with physical layer cell identity NCellID, when
        % NRSID field is not present or is empty
        dmrs.NRSID = carrier.NCellID;
    end
    % Get the value of NumCDMGroupsWithoutData property depending on the
    % value of field NumCDMGroupsWithoutData
    if isfield(chs,'NumCDMGroupsWithoutData') && ~isempty(chs.NumCDMGroupsWithoutData) ...
            && ~pusch.TransformPrecoding && chs.NumCDMGroupsWithoutData
        dmrs.NumCDMGroupsWithoutData = chs.NumCDMGroupsWithoutData;
    else
        % When the field NumCDMGroupsWithoutData is not present or is set
        % to empty or zero, assign the value to 1 or 2, depending on
        % transform precoding disabled or enabled
        dmrs.NumCDMGroupsWithoutData = 1 + pusch.TransformPrecoding;
    end
    pusch.DMRS = dmrs;

    % Set PT-RS parameters
    pusch = passign(chs,pusch,'EnablePTRS');
    ptrs = nrPUSCHPTRSConfig;
    ptrs = passign(chs,ptrs,'PTRSTimeDensity','TimeDensity');
    ptrs = passign(chs,ptrs,'PTRSFrequencyDensity','FrequencyDensity',~pusch.TransformPrecoding);
    ptrs = passign(chs,ptrs,'PTRSNumSamples','NumPTRSSamples',pusch.TransformPrecoding);
    ptrs = passign(chs,ptrs,'PTRSNumGroups','NumPTRSGroups',pusch.TransformPrecoding);
    ptrs = passign(chs,ptrs,'PTRSREOffset','REOffset');
    ptrs = passign(chs,ptrs,'PTRSPortSet');
    [ptrs,fieldPresent] = passign(chs,ptrs,'PTRSNID','NID',pusch.TransformPrecoding);
    if ~fieldPresent
        % Assign NID with DM-RS scrambling identity NRSID, when PTRSNID
        % field is not present or is empty
        ptrs.NID = pusch.DMRS.NRSID;
    end
    pusch.PTRS = ptrs;

    % Disable PT-RS
    if pusch.TransformPrecoding
        cond1 = isfield(chs,'PTRSNumSamples') && isempty(chs.PTRSNumSamples);
        cond2 = isfield(chs,'PTRSNumGroups') && isempty(chs.PTRSNumGroups);
        cond = cond1 | cond2;
    else
        cond = isfield(chs,'PTRSFrequencyDensity') && isempty(chs.PTRSFrequencyDensity);
    end
    if cond || (isfield(chs,'PTRSTimeDensity') && isempty(chs.PTRSTimeDensity))
        pusch.EnablePTRS = 0;
    end

end

function [o,cond] = passign(s,o,f,p,ac)

    cond = isfield(s,f) && ~isempty(s.(f));
    if nargin == 5
        cond = cond && ac;
    end

    if cond
        if nargin == 3
            o.(f) = s.(f);
        else
            o.(p) = s.(f);
        end
    end
end
