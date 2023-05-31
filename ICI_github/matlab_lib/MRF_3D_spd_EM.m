function [PI_IN, S_POST, beta_upd, alpha_upd] = MRF_3D_spd_EM(obj, PI_OUT)
            
            % Get sizing information
            N = size(PI_OUT, 1);
            
            % Enclose the 3-space defined by the coordinates property
            % within the smallest containing cube
            MaxVals = max(obj.coordinates, [], 1);
            Nx = MaxVals(1); Ny = MaxVals(2); Nz = MaxVals(3);
            
            
            % Declare model parameter-related constants        
            ebx = exp(obj.betax);
            eby = exp(obj.betay);
            ebz = exp(obj.betaz);
            ebxi = exp(-obj.betax);
            ebyi = exp(-obj.betay);
            ebzi = exp(-obj.betaz);
            e0 = exp(obj.alpha);
            e1 = exp(-obj.alpha);
            
            % checkerboard pattern indices in 3-space
            chk = NaN(Nx,Ny,Nz);
            check = checkerboard(1, Nx, Ny);
            check = check(1:Nx,1:Ny) > 0;
            for i = 1:Nz
                if mod(i, 2) == 0
                    chk(:,:,i) = ~check;
                else
                    chk(:,:,i) = check;
                end
            end
            blackIdx = find(chk == 0);
            whiteIdx = find(chk == 1);
            % Note the linear indices that belong to actual signal
            % coordinates...
            presentIdx = sub2ind([Nx, Ny, Nz], obj.coordinates(:,1), ...
                obj.coordinates(:,2), obj.coordinates(:,3));
            % ...and those that belong to dummy nodes
            missingIdx = setdiff(1:Nx*Ny*Nz, presentIdx);
            % Keep only those indices that were assigned coordinates
            blackIdx = intersect(blackIdx, presentIdx);
            whiteIdx = intersect(whiteIdx, presentIdx);
            
            
            % Initialize incoming probability messages
            nodePot = 1/2 * ones(Nx,Ny,Nz,2);
            nodePot(presentIdx) = 1 - PI_OUT;       	% Prob. that node is 0
            nodePot(presentIdx + Nx*Ny*Nz) = PI_OUT;	% Prob. that node is 1
            
            % Initialize messages
            msgFromRight = 1/2 * ones(Nx,Ny,Nz,2);
            msgFromLeft = 1/2 * ones(Nx,Ny,Nz,2);
            msgFromTop = 1/2 * ones(Nx,Ny,Nz,2);
            msgFromBottom = 1/2 * ones(Nx,Ny,Nz,2);
            msgFromFront = 1/2 * ones(Nx,Ny,Nz,2);
            msgFromBack = 1/2 * ones(Nx,Ny,Nz,2);
            
            prod0 = zeros(Nx,Ny,Nz);
            prod1 = prod0;
            
            for iter = 1:2*obj.maxIter
                % First grab indices of messages to update this round
                if(mod(iter,2) == 1)
                    ind = blackIdx;
                    ind1 = blackIdx;
                    ind2 = blackIdx + Nx*Ny*Nz;     % Offset linear idx
                else
                    ind = whiteIdx;
                    ind1 = whiteIdx;
                    ind2 = whiteIdx + Nx*Ny*Nz;     % Offset linear idx
                end
                
                % Update messages from left 
                prod0(:,2:end,:) = e0*nodePot(:,1:end-1,:,1) .* ...
                    msgFromLeft(:,1:end-1,:,1) .* msgFromTop(:,1:end-1,:,1) .* ...
                    msgFromBottom(:,1:end-1,:,1) .* msgFromFront(:,1:end-1,:,1) .* ...
                    msgFromBack(:,1:end-1,:,1);
                prod1(:,2:end,:) = e1*nodePot(:,1:end-1,:,2) .* ...
                    msgFromLeft(:,1:end-1,:,2) .* msgFromTop(:,1:end-1,:,2) .* ...
                    msgFromBottom(:,1:end-1,:,2) .* msgFromFront(:,1:end-1,:,2) .* ...
                    msgFromBack(:,1:end-1,:,2);
                %%
                prod0(:,1,:) = e0*nodePot(:,end,:,1) .* ...
                    msgFromLeft(:,end,:,1) .* msgFromTop(:,end,:,1) .* ...
                    msgFromBottom(:,end,:,1) .* msgFromFront(:,end,:,1) .* ...
                    msgFromBack(:,end,:,1);
                prod1(:,1,:) = e1*nodePot(:,end,:,2) .* ...
                    msgFromLeft(:,end,:,2) .* msgFromTop(:,end,:,2) .* ...
                    msgFromBottom(:,end,:,2) .* msgFromFront(:,end,:,2) .* ...
                    msgFromBack(:,end,:,2);
                %%
                p0 = prod0*ebx + prod1*ebxi;
                p1 = prod0*ebxi + prod1*ebx;
                sump0p1 = p0+p1;
                
                msgFromLeft(ind1) = p0(ind) ./ sump0p1(ind);
                msgFromLeft(ind2) = p1(ind) ./ sump0p1(ind);
%                 msgFromLeft(:,1,:,:) = 1/2;     	% Dummy edge msgs
                % The messages that are arriving from the left at nodes
                % whose left-side neighbors are dummy nodes should be
                % uninformative (1/2), thus we must correct such messages.
                % We can identify these nodes that have left-side dummy
                % neighbors using some linear indexing tricks
%                 correctedIdx = missingIdx + Nx;
%                 correctedIdx = correctedIdx(correctedIdx <= Nx*Ny*Nz);
%                 correctedIdx = [correctedIdx; correctedIdx + Nx*Ny*Nz];
%                 msgFromLeft(correctedIdx) = 1/2;    % Dummy missing msgs
                
                
                 % Update messages from right 
                prod0(:,1:end-1,:) = e0*nodePot(:,2:end,:,1) .* ...
                    msgFromRight(:,2:end,:,1) .* msgFromTop(:,2:end,:,1) .* ...
                    msgFromBottom(:,2:end,:,1) .* msgFromFront(:,2:end,:,1) .* ...
                    msgFromBack(:,2:end,:,1);
                prod1(:,1:end-1,:) = e1*nodePot(:,2:end,:,2) .* ...
                    msgFromRight(:,2:end,:,2) .* msgFromTop(:,2:end,:,2) .* ...
                    msgFromBottom(:,2:end,:,2) .* msgFromFront(:,2:end,:,2) .* ...
                    msgFromBack(:,2:end,:,2);
                %%
                prod0(:,end,:) = e0*nodePot(:,1,:,1) .* ...
                    msgFromRight(:,1,:,1) .* msgFromTop(:,1,:,1) .* ...
                    msgFromBottom(:,1,:,1) .* msgFromFront(:,1,:,1) .* ...
                    msgFromBack(:,1,:,1);
                prod1(:,end,:) = e1*nodePot(:,1,:,2) .* ...
                    msgFromRight(:,1,:,2) .* msgFromTop(:,1,:,2) .* ...
                    msgFromBottom(:,1,:,2) .* msgFromFront(:,1,:,2) .* ...
                    msgFromBack(:,1,:,2);
                p0 = prod0*ebx + prod1*ebxi;
                p1 = prod0*ebxi + prod1*ebx;
                sump0p1 = p0 + p1;
                
                msgFromRight(ind1) = p0(ind) ./ sump0p1(ind);
                msgFromRight(ind2) = p1(ind) ./ sump0p1(ind);
%                 msgFromRight(:,end,:,:) = 1/2;    	% Dummy edge msgs
                % As in the previous case, we must manually correct
                % messages at nodes whose right-side neighbors are dummy
                % nodes. Again, use linear indexing tricks
%                 correctedIdx = missingIdx - Nx;
%                 correctedIdx = correctedIdx(correctedIdx >= 1);
%                 correctedIdx = [correctedIdx; correctedIdx + Nx*Ny*Nz];
%                 msgFromRight(correctedIdx) = 1/2;   % Dummy missing msgs
                
                
                % Update messages from top 
                prod0(2:end,:,:) = e0*nodePot(1:end-1,:,:,1) .* ...
                    msgFromLeft(1:end-1,:,:,1) .* msgFromTop(1:end-1,:,:,1) .* ...
                    msgFromRight(1:end-1,:,:,1) .* msgFromFront(1:end-1,:,:,1) .* ...
                    msgFromBack(1:end-1,:,:,1);
                prod1(2:end,:,:) = e1*nodePot(1:end-1,:,:,2) .* ...
                    msgFromLeft(1:end-1,:,:,2) .* msgFromTop(1:end-1,:,:,2) .* ...
                    msgFromRight(1:end-1,:,:,2) .* msgFromFront(1:end-1,:,:,2) .* ...
                    msgFromBack(1:end-1,:,:,2);
                %%
                prod0(1,:,:) = e0*nodePot(end,:,:,1) .* ...
                    msgFromLeft(end,:,:,1) .* msgFromTop(end,:,:,1) .* ...
                    msgFromRight(end,:,:,1) .* msgFromFront(end,:,:,1) .* ...
                    msgFromBack(end,:,:,1);
                prod1(1,:,:) = e1*nodePot(end,:,:,2) .* ...
                    msgFromLeft(end,:,:,2) .* msgFromTop(end,:,:,2) .* ...
                    msgFromRight(end,:,:,2) .* msgFromFront(end,:,:,2) .* ...
                    msgFromBack(end,:,:,2);
                p0 = prod0*eby + prod1*ebyi;
                p1 = prod0*ebyi + prod1*eby;
                sump0p1 = p0 + p1;
                
                msgFromTop(ind1) = p0(ind) ./ sump0p1(ind);
                msgFromTop(ind2) = p1(ind) ./ sump0p1(ind);
%                 msgFromTop(1,:,:,:) = 1/2;          % Dummy edge msgs
                % As in the previous case, we must manually correct
                % messages at nodes whose top-side neighbors are dummy
                % nodes. Again, use linear indexing tricks
%                 correctedIdx = missingIdx + 1;
%                 correctedIdx = correctedIdx(correctedIdx <= Nx*Ny*Nz);
%                 correctedIdx = [correctedIdx; correctedIdx + Nx*Ny*Nz];
%                 msgFromTop(correctedIdx) = 1/2;  	% Dummy missing msgs
                
                
                % Update messages from bottom 
                prod0(1:end-1,:,:) = e0*nodePot(2:end,:,:,1) .* ...
                    msgFromRight(2:end,:,:,1) .* msgFromLeft(2:end,:,:,1) .* ...
                    msgFromBottom(2:end,:,:,1) .* msgFromFront(2:end,:,:,1) .* ...
                    msgFromBack(2:end,:,:,1);
                prod1(1:end-1,:,:) = e1*nodePot(2:end,:,:,2) .* ...
                    msgFromRight(2:end,:,:,2) .* msgFromLeft(2:end,:,:,2) .* ...
                    msgFromBottom(2:end,:,:,2) .* msgFromFront(2:end,:,:,2) .* ...
                    msgFromBack(2:end,:,:,2);
                %%
                prod0(end,:,:) = e0*nodePot(1,:,:,1) .* ...
                    msgFromRight(1,:,:,1) .* msgFromLeft(1,:,:,1) .* ...
                    msgFromBottom(1,:,:,1) .* msgFromFront(1,:,:,1) .* ...
                    msgFromBack(1,:,:,1);
                prod1(end,:,:) = e1*nodePot(1,:,:,2) .* ...
                    msgFromRight(1,:,:,2) .* msgFromLeft(1,:,:,2) .* ...
                    msgFromBottom(1,:,:,2) .* msgFromFront(1,:,:,2) .* ...
                    msgFromBack(1,:,:,2);
                p0 = prod0*eby + prod1*ebyi;
                p1 = prod0*ebyi + prod1*eby;
                sump0p1 = p0 + p1;
                
                msgFromBottom(ind1) = p0(ind) ./ sump0p1(ind);
                msgFromBottom(ind2) = p1(ind) ./ sump0p1(ind);
%                 msgFromBottom(end,:,:,:) = 1/2;     % Dummy edge msgs
                % As in the previous case, we must manually correct
                % messages at nodes whose bottom-side neighbors are dummy
                % nodes. Again, use linear indexing tricks
%                 correctedIdx = missingIdx - 1;
%                 correctedIdx = correctedIdx(correctedIdx >= 1);
%                 correctedIdx = [correctedIdx; correctedIdx + Nx*Ny*Nz];
%                 msgFromBottom(correctedIdx) = 1/2;	% Dummy missing msgs
                
                
                % Update messages from front 
                prod0(:,:,2:end) = e0*nodePot(:,:,1:end-1,1) .* ...
                    msgFromLeft(:,:,1:end-1,1) .* msgFromTop(:,:,1:end-1,1) .* ...
                    msgFromRight(:,:,1:end-1,1) .* msgFromFront(:,:,1:end-1,1) .* ...
                    msgFromBottom(:,:,1:end-1,1);
                prod1(:,:,2:end) = e1*nodePot(:,:,1:end-1,2) .* ...
                    msgFromLeft(:,:,1:end-1,2) .* msgFromTop(:,:,1:end-1,2) .* ...
                    msgFromRight(:,:,1:end-1,2) .* msgFromFront(:,:,1:end-1,2) .* ...
                    msgFromBottom(:,:,1:end-1,2);
                %%
                prod0(:,:,1) = e0*nodePot(:,:,end,1) .* ...
                    msgFromLeft(:,:,end,1) .* msgFromTop(:,:,end,1) .* ...
                    msgFromRight(:,:,end,1) .* msgFromFront(:,:,end,1) .* ...
                    msgFromBottom(:,:,end,1);
                prod1(:,:,1) = e1*nodePot(:,:,end,2) .* ...
                    msgFromLeft(:,:,end,2) .* msgFromTop(:,:,end,2) .* ...
                    msgFromRight(:,:,end,2) .* msgFromFront(:,:,end,2) .* ...
                    msgFromBottom(:,:,end,2);
                p0 = prod0*ebz + prod1*ebzi;
                p1 = prod0*ebzi + prod1*ebz;
                sump0p1 = p0 + p1;
                
                msgFromFront(ind1) = p0(ind) ./ sump0p1(ind);
                msgFromFront(ind2) = p1(ind) ./ sump0p1(ind);
%                 msgFromFront(:,:,1,:) = 1/2;        % Dummy edge msgs
                % As in the previous case, we must manually correct
                % messages at nodes whose front-side neighbors are dummy
                % nodes. Again, use linear indexing tricks
%                 correctedIdx = missingIdx + Nx*Ny;
%                 correctedIdx = correctedIdx(correctedIdx <= Nx*Ny*Nz);
%                 correctedIdx = [correctedIdx; correctedIdx + Nx*Ny*Nz];
%                 msgFromFront(correctedIdx) = 1/2;	% Dummy missing msgs
                
                
                % Update messages from back 
                prod0(:,:,1:end-1) = e0*nodePot(:,:,2:end,1) .* ...
                    msgFromRight(:,:,2:end,1) .* msgFromLeft(:,:,2:end,1) .* ...
                    msgFromBottom(:,:,2:end,1) .* msgFromTop(:,:,2:end,1) .* ...
                    msgFromBack(:,:,2:end,1);
                prod1(:,:,1:end-1) = e1*nodePot(:,:,2:end,2) .* ...
                    msgFromRight(:,:,2:end,2) .* msgFromLeft(:,:,2:end,2) .* ...
                    msgFromBottom(:,:,2:end,2) .* msgFromTop(:,:,2:end,2) .* ...
                    msgFromBack(:,:,2:end,2);
                %%
                prod0(:,:,end) = e0*nodePot(:,:,1,1) .* ...
                    msgFromRight(:,:,1,1) .* msgFromLeft(:,:,1,1) .* ...
                    msgFromBottom(:,:,1,1) .* msgFromTop(:,:,1,1) .* ...
                    msgFromBack(:,:,1,1);
                prod1(:,:,end) = e1*nodePot(:,:,1,2) .* ...
                    msgFromRight(:,:,1,2) .* msgFromLeft(:,:,1,2) .* ...
                    msgFromBottom(:,:,1,2) .* msgFromTop(:,:,1,2) .* ...
                    msgFromBack(:,:,1,2);
                p0 = prod0*ebz + prod1*ebzi;
                p1 = prod0*ebzi + prod1*ebz;
                sump0p1 = p0 + p1;
                
                msgFromBack(ind1) = p0(ind) ./ sump0p1(ind);
                msgFromBack(ind2) = p1(ind) ./ sump0p1(ind);
%                 msgFromBack(:,:,end,:) = 1/2;       % Dummy edge msgs
                % As in the previous case, we must manually correct
                % messages at nodes whose back-side neighbors are dummy
                % nodes. Again, use linear indexing tricks
%                 correctedIdx = missingIdx - Nx*Ny;
%                 correctedIdx = correctedIdx(correctedIdx >= 1);
%                 correctedIdx = [correctedIdx; correctedIdx + Nx*Ny*Nz];
%                 msgFromBack(correctedIdx) = 1/2;    % Dummy missing msgs
                
            end
            
            
            % Compute extrinsic likelihood, marginal potential and s_hat
            msgProds = msgFromLeft .* msgFromRight .* msgFromTop .* ...
                msgFromBottom .* msgFromFront .* msgFromBack;
            msgProds(:,:,:,1) = msgProds(:,:,:,1)*e0;
            msgProds(:,:,:,2) = msgProds(:,:,:,2)*e1;
            Le_spdec = log(msgProds(:,:,:,2)./msgProds(:,:,:,1));
            PI_IN = 1 ./ (1 + exp(-Le_spdec(presentIdx)));
            
            msgProds = msgProds.*nodePot;
            sumMsgProds = sum(msgProds, 4);
            S_POST = msgProds(:,:,:,2) ./ sumMsgProds;    % Pr{S(n) = 1 | Y}
            S_POST = S_POST(presentIdx);
            S_HAT = double(S_POST > 1/2);
            
            PI = msgProds(:,:,:,2) ./ sumMsgProds;    % Pr{S(n) = 1 | Y}
            
            
            
            
            % Compute parameter updates (will require MATLAB's Optimization
            % Toolbox).  Currently learns a single value of beta
            beta_upd = [obj.betax, obj.betay, obj.betaz];
            alpha_upd = obj.alpha;
            options = optimset('GradObj', 'on', 'Hessian', ...
                'off', 'MaxFunEvals', 1000, 'tolX', 1e-14, 'Display', ...
                'notify', 'Algorithm', 'interior-point');
            lb = [0; 0; 0; 0];   % Lower bounds [alpha; beta]
            ub = [10; 30; 30; 30];    % Upper bounds [alpha; beta]
            if obj.EMcnt >= 1 
                if strcmpi(obj.learn_alpha, 'true') || ...
                        strcmpi(obj.learn_beta, 'true')
                    [updates, ~, exitFlag] = fmincon(@meanfieldLF_3D, ...
                        [obj.alpha; obj.betax;obj.betay;obj.betaz], [], [], [], [], lb, ub, [], ...
                        options, PI);
    %                 fprintf('exitflag=%d\n',exitFlag);
                    if (exitFlag>0)
                        alpha_upd = updates(1);
                        beta_upd = updates(2:end);
                    end
                end
%                 if strcmpi(obj.learn_alpha, 'true'), alpha_upd = updates(1); end
%                 if strcmpi(obj.learn_beta, 'true'), beta_upd = updates(2:end); end
            end
            obj.EMcnt = obj.EMcnt + 1;

%             % *********************************
%             % No EM updates for the time being
%             % *********************************
%             beta_upd = mean([obj.betax, obj.betay, obj.betaz]);
%             alpha_upd = obj.alpha; 
            
end