function [pathGains,pathFilters,offset,sampleTimes,ofdminfo,sampleTimesDefaulted,initialNSlot,hasSampleRate] = my_ofdm_info(varargin)
    
    coder.extrinsic('nr5g.internal.parseOptions');
    fcnName = 'nrPerfectChannelEstimate';
    
    % Determine if syntax with nrCarrierConfig is being used and parse
    % relevant inputs
    isCarrierSyntax = isa(varargin{1},'nrCarrierConfig');
    if (isCarrierSyntax) % CARRIER,PATHGAINS,PATHFILTERS,...
        carrier = varargin{1};
        validateattributes(carrier,{'nrCarrierConfig'}, ...
            {'scalar'},fcnName,'Carrier specific configuration object');
        pathGains = varargin{2};
        pathFilters = varargin{3};
        initialNSlot = carrier.NSlot;
        optstart = 4;
    else % PATHGAINS,PATHFILTERS,NRB,SCS,INITIALNSLOT,...
        narginchk(5,17);
        pathGains = varargin{1};
        pathFilters = varargin{2};
        NRB = varargin{3};
        SCS = varargin{4};
        initialNSlot = varargin{5};
        optstart = 6;
    end
    
    % Validate channel path gains
    validateattributes(pathGains,{'double','single'}, ...
        {},fcnName,'PATHGAINS');
    coder.internal.errorIf(ndims(pathGains)>4, ...
        'nr5g:nrPerfectChannelEstimate:InvalidPathDims',ndims(pathGains));
    
    % Validate path filters impulse response
    validateattributes(pathFilters,{'double'}, ...
        {'2d'},fcnName,'PATHFILTERS');
    coder.internal.errorIf(size(pathGains,2)~=size(pathFilters,2), ...
        'nr5g:nrPerfectChannelEstimate:InconsistentPaths', ...
        size(pathGains,2),size(pathFilters,2));
    
    if (~isCarrierSyntax)
        % Validate the number of resource blocks (1...275)
        validateattributes(NRB,{'numeric'}, ...
            {'real','integer','scalar','>=',1,'<=',275},fcnName,'NRB');

        % Validate subcarrier spacing input in kHz (15/30/60/120/240)
        validateattributes(SCS,{'numeric'}, ...
            {'real','integer','scalar'},fcnName,'SCS');
        validSCS = [15 30 60 120 240];
        coder.internal.errorIf(~any(SCS==validSCS), ...
            'nr5g:nrPerfectChannelEstimate:InvalidSCS', ...
            SCS,num2str(validSCS));

        % Validate zero-based initial slot number
        validateattributes(initialNSlot,{'numeric'}, ...
            {'real','nonnegative','scalar','integer'}, ...
            fcnName,'INITIALNSLOT');
    end
    
    % Parse optional arguments: 'offset', 'sampleTimes', and options.
    % Options are either name-value pairs, or a char/string value for the
    % cyclic prefix length. These options can appear in any position among
    % the 'offset' and 'sampleTimes' arguments, but multiple name-value
    % pairs must appear together. The value-only cyclic prefix length is
    % only valid for syntaxes without nrCarrierConfig
    argpos = [0 0]; % positions of offset and sampleTimes, 0 if absent
    pos = 1;
    firstoptarg = 0; % position of first option, 0 if absent
    cpnamefn = @(x)any(strcmpi(x,{'normal','extended'}));
    for i = optstart:nargin
        if (ischar(varargin{i}) || isstring(varargin{i}) ...
                || isstruct(varargin{i}) || pos>2)
            if (firstoptarg==0)
                firstoptarg = i;
            end
        else
            if (firstoptarg==0 || cpnamefn(varargin{firstoptarg}) ...
                    || isstruct(varargin{firstoptarg}) || ...
                    (mod(i,2) == mod(firstoptarg,2)))
                argpos(pos) = i;
                pos = pos + 1;
            end
        end
    end
    
    if (firstoptarg~=0)
        lastoptarg = min([argpos(argpos~=0 & argpos>firstoptarg)-1 nargin]);
    end
    
    % Parse options and get OFDM information
    if (isCarrierSyntax)
        optNames = {'Nfft','SampleRate','CyclicPrefixFraction'};
        if firstoptarg~=0
            opts = coder.const(nr5g.internal.parseOptions([fcnName '(carrier,...'],optNames,varargin{firstoptarg:lastoptarg}));
        else
            optins = {};
            opts = coder.const(nr5g.internal.parseOptions([fcnName '(carrier,...'],optNames,optins{:}));
        end
        ofdminfo = nr5g.internal.OFDMInfo(carrier,opts);
    else
        optNames = {'CyclicPrefix','Nfft','SampleRate','CyclicPrefixFraction'};
        if firstoptarg~=0
            if (numel(varargin)>firstoptarg && ~cpnamefn(varargin{firstoptarg}))
                opts = coder.const(nr5g.internal.parseOptions(fcnName,optNames,varargin{firstoptarg:lastoptarg}));
                ECP = strcmpi(opts.CyclicPrefix,'extended');
                % If there are NV pairs, ofdminfo must be constant
                ofdminfo = coder.const(feval('nr5g.internal.OFDMInfo',NRB,SCS,ECP,opts));
            else
                coder.internal.errorIf(lastoptarg~=firstoptarg, ...
                    'nr5g:nrPerfectChannelEstimate:CPValueOnlyOption');
                cp = varargin{firstoptarg};
                validateattributes(cp,{'char' 'string'},{},fcnName,'CP');
                cp = validatestring(cp,{'normal','extended'},fcnName,'CP');
                ECP = strcmpi(cp,'extended');
                % This allows changing the cp at runtime
                opts = coder.const(nr5g.internal.parseOptions(fcnName,optNames));
                if ECP
                    opts.CyclicPrefix = 'extended';
                else
                    opts.CyclicPrefix = 'normal';
                end
                ofdminfo = nr5g.internal.OFDMInfo(NRB,SCS,ECP,opts);
            end
        else
            opts = coder.const(nr5g.internal.parseOptions(fcnName,optNames));
            ofdminfo = nr5g.internal.OFDMInfo(NRB,SCS,0,opts);
        end
    end
    
    % If performing code generation, the presence of sample rate with the
    % function syntax using nrCarrierConfig triggers a compile-time error
    hasSampleRate = ~isempty(opts.SampleRate);
    if (~isempty(coder.target) && isCarrierSyntax && hasSampleRate)
        coder.internal.errorIf(true,'nr5g:nrPerfectChannelEstimate:CompilationCarrierSampleRate',sprintf('%g',opts.SampleRate),'IfNotConst','Fail');
    end
    
    % Validate offset or provide a default value
    if argpos(1)~=0
        offset = varargin{argpos(1)};
        Nh = size(pathFilters,1);
        validateattributes(offset,{'numeric'}, ...
            {'real','nonnegative','scalar'},fcnName,'offset');
        coder.internal.errorIf(offset>(Nh-1), ...
            'nr5g:nrPerfectChannelEstimate:InvalidOffset', ...
            offset,Nh);
        mincplen = min(ofdminfo.CyclicPrefixLengths);
        coder.internal.errorIf(offset>mincplen, ...
            'nr5g:nrPerfectChannelEstimate:InvalidOffsetCPLength', ...
            offset,mincplen);
    else
        % Default: use nrPerfectTimingEstimate to establish the timing
        % offset
        offset = nrPerfectTimingEstimate(pathGains,pathFilters);
    end
    
    % Validate sampleTimes or provide a default value
    if argpos(2)~=0
        sampleTimes = varargin{argpos(2)};
        Ncs = size(pathGains,1);
        validateattributes(sampleTimes,{'double'}, ...
            {'column','increasing'},fcnName,'sampleTimes');
        coder.internal.errorIf(length(sampleTimes)~=Ncs, ...
            'nr5g:nrPerfectChannelEstimate:InvalidSampleTimes', ...
            length(sampleTimes),Ncs);
        sampleTimesDefaulted = false;
    else
        % Default: vector of times at the OFDM sampling rate, one for each
        % channel snapshot in 'pathGains' and starting at zero
        sampleTimes = (0:(size(pathGains,1)-1)).' / ofdminfo.SampleRate;
        sampleTimesDefaulted = true;
    end
    
    ofdminfo.NSubcarriers=ofdminfo.Nfft;
    
end