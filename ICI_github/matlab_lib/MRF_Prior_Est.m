function [pi_out, S_POST] = MRF_Prior_Est(pi_in_alpha,A,B, Iteration_For_MRF)

[N,M] = size(pi_in_alpha);
nodePotFunction = zeros(N,M,2);

nodePotFunction(:,:,1) = 1 - pi_in_alpha;    % Prob. that node is 0
nodePotFunction(:,:,2) = pi_in_alpha;        % Prob. that node is 1

% Initialize messages, msg(:,:,1) denote message for alpha = -1, and
% msg(:,:,2) denote message for alpha = 1;
msgFromRight = 0.5*ones(N,M,2);
msgFromLeft = 0.5*ones(N,M,2);
msgFromTop = 0.5*ones(N,M,2);
msgFromBottom = 0.5*ones(N,M,2);

% message matrix initialize. 
prob0 = zeros(N,M);
prob1 = zeros(N,M);

% blink schedule, update nodes iteratively:
checker_mat = checkerboard(1, N/2, M/2);
K1 = find(reshape(checker_mat,[],1) == 0);
K2 = find(reshape(checker_mat,[],1) > 0);
blink_vec = zeros(N*M,1);
blink_vec(K2) = 1;
blink_vec(K1) = 0;
blink_mat = reshape(blink_vec,N,M);
[upt_1_x, upt_1_y] = find(blink_mat == 0);
[upt_2_x, upt_2_y] = find(blink_mat == 1);

for iter = 1:Iteration_For_MRF
    
    if(mod(iter,2) == 1)
        x = upt_1_x; y = upt_1_y;
    else
        x = upt_2_x; y = upt_2_y;
    end
    idx = sub2ind([N, M], x, y);
    idx1 = sub2ind([N, M, 2], x, y, ones(numel(x),1)); % correspond to (:,:,1)
    idx2 = sub2ind([N, M, 2], x, y, 2*ones(numel(x),1)); % correspond to (:,:,2)
    
    % update messages from left
    % 1. message for alpha = -1
    prob0(:,2:end) = exp(A) * (1-pi_in_alpha(:,1:end-1)) .* ...
                    msgFromLeft(:,1:end-1,1) .* msgFromTop(:,1:end-1,1) .* msgFromBottom(:,1:end-1,1);
                
    % 2. message for alpha = 1
    prob1(:,2:end) = exp(-A) * (pi_in_alpha(:,1:end-1)) .* ...
                    msgFromLeft(:,1:end-1,2) .* msgFromTop(:,1:end-1,2) .* msgFromBottom(:,1:end-1,2);
    % 3. compute probability of alpha = -1/1
    p_0 = prob0*exp(B) + prob1*exp(-B);
    p_1 = prob0*exp(-B) + prob1*exp(B);           
    % 4. normalized them
    p0 = p_0(idx)./(p_0(idx) + p_1(idx));
    p1 = p_1(idx)./(p_0(idx) + p_1(idx));
    
    msgFromLeft(idx1) = p0;
    msgFromLeft(idx2) = p1;
    msgFromLeft(:,1,:) = 0.5;
    
    % update messages from right
    % 1. message for alpha = -1
    prob0(:,1:end-1) = exp(A) * (1-pi_in_alpha(:,2:end)) .* ...
                    msgFromRight(:,2:end,1) .* msgFromTop(:,2:end,1) .* msgFromBottom(:,2:end,1);
                
    % 2. message for alpha = 1
    prob1(:,1:end-1) = exp(-A) * (pi_in_alpha(:,2:end)) .* ...
                    msgFromRight(:,2:end,2) .* msgFromTop(:,2:end,2) .* msgFromBottom(:,2:end,2);
    % 3. compute probability of alpha = -1/1
    p_0 = prob0*exp(B) + prob1*exp(-B);
    p_1 = prob0*exp(-B) + prob1*exp(B);           
    % 4. normalized them
    p0 = p_0(idx)./(p_0(idx) + p_1(idx));
    p1 = p_1(idx)./(p_0(idx) + p_1(idx));
    
    msgFromRight(idx1) = p0;
    msgFromRight(idx2) = p1;
    msgFromRight(:,end,:) = 0.5;
    
    
    % update messages from top
    % 1. message for alpha = -1
    prob0(2:end,:) = exp(A) * (1-pi_in_alpha(1:end-1,:)) .* ...
                    msgFromLeft(1:end-1,:,1) .* msgFromTop(1:end-1,:,1) .* msgFromRight(1:end-1,:,1);
                
    % 2. message for alpha = 1
    prob1(2:end,:) = exp(-A) * (pi_in_alpha(1:end-1,:)) .* ...
                    msgFromLeft(1:end-1,:,2) .* msgFromTop(1:end-1,:,2) .* msgFromRight(1:end-1,:,2);
    % 3. compute probability of alpha = -1/1
    p_0 = prob0*exp(B) + prob1*exp(-B);
    p_1 = prob0*exp(-B) + prob1*exp(B);           
    % 4. normalized them
    p0 = p_0(idx)./(p_0(idx) + p_1(idx));
    p1 = p_1(idx)./(p_0(idx) + p_1(idx));
    
    msgFromTop(idx1) = p0;
    msgFromTop(idx2) = p1;
    msgFromTop(1,:,:) = 0.5;
    
    % update messages from bottom
    % 1. message for alpha = -1
    prob0(1:end-1,:) = exp(A) * (1-pi_in_alpha(2:end,:)) .* ...
                    msgFromLeft(2:end,:,1) .* msgFromBottom(2:end,:,1) .* msgFromRight(2:end,:,1);
                
    % 2. message for alpha = 1
    prob1(1:end-1,:) = exp(-A) * (pi_in_alpha(2:end,:)) .* ...
                    msgFromLeft(2:end,:,2) .* msgFromBottom(2:end,:,2) .* msgFromRight(2:end,:,2);
    % 3. compute probability of alpha = -1/1
    p_0 = prob0*exp(B) + prob1*exp(-B);
    p_1 = prob0*exp(-B) + prob1*exp(B);           
    % 4. normalized them
    p0 = p_0./(p_0 + p_1);
    p1 = p_1./(p_0 + p_1);
    
    p0 = p_0(idx)./(p_0(idx) + p_1(idx));
    p1 = p_1(idx)./(p_0(idx) + p_1(idx));
    
    msgFromBottom(idx1) = p0;
    msgFromBottom(idx2) = p1;
    msgFromBottom(end,:,:) = 0.5;
    
end

msg_gather = msgFromLeft .* msgFromRight .* msgFromTop .* msgFromBottom;
msgProds(:,:,1) = msg_gather(:,:,1)*exp(A);
msgProds(:,:,2) = msg_gather(:,:,2)*exp(-A);

pi_out = squeeze(msgProds(:,:,2) ./ (msgProds(:,:,1)+msgProds(:,:,2)));


NodeGatherMsg = msgProds.*nodePotFunction;
sumNodeGatherMsg = sum(NodeGatherMsg,3);

S_POST = NodeGatherMsg(:,:,2) ./ sumNodeGatherMsg;    % Pr{S(n,t) = 1 | Y}

% S_HAT = double(S_POST > 1/2);
end

