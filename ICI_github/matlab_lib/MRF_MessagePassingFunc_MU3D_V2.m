function [NMSE,NMSE_individual,X_Est,DebugData,sysParComm] = MRF_MessagePassingFunc_MU3D_V2(X_org, A, Y, sysPar, sysParComm, Iter_max, debug_mode, Prior_Model)
% X_org: NxM size original channel matrix
% A: PxNxK Sensing Matrix for different users. 
% Y: PxMxK received signal for K users.  
% sysPar: structure of system parameters
% Damling adoptted, 3D Ising model adoptted, Hidden-MRF is also optional. 
% In this version, the LMMSE is remade to be more general with snesing
% matrix
MRF_Inner_Iteration = 40;
[N,M,K] = size(X_org);
% S = N*M;
[Tp,~,~] = size(Y); %
% binary MRF prior parameters, modified later
a = sysParComm.a;% larger a encourage sparser matrix
b = sysParComm.b;% larger b encourage larger non-zero clusters. 
c = sysParComm.c;% larger c encourage larger non-zero clusters. 
beta = sysParComm.damping_par;
EM_learning = sysParComm.EM_enable;
X_A_PRI = zeros(N,M,K);
parfor user_cnt = 1:K
    V_A_PRI(:,user_cnt) = sysPar(user_cnt).sigVar;
end
X_B_PRI = zeros(N,M,K);
V_B_PRI = zeros(M,K);
X_A_POST = zeros(N,M,K);
V_A_POST = zeros(M,K);
X_B_POST = zeros(N,M,K);
V_B_POST = zeros(M,K);

X_A_EXT = zeros(N,M,K);
V_A_EXT = zeros(M,K);
X_B_EXT = zeros(N,M,K);
V_B_EXT = zeros(M,K);

X_A_EXT_Last = zeros(N,M,K);
V_A_EXT_Last = zeros(M,K);
X_B_EXT_Last = zeros(N,M,K);
V_B_EXT_Last = zeros(M,K);

X_Est = zeros(N,M,K);

pi_in_alpha = zeros(N,M,K);
pi_out_alpha = zeros(N,M,K);


for iter_cnt = 1:Iter_max
    
    for user_cnt = 1:K
        
%         Y_k = squeeze(Y(:,:,user_cnt));
        Y_k = Y;
        A_k = squeeze(A(:,:,user_cnt));
        sigma = sysPar(user_cnt).noiseVar;
        
        for ant_idx = 1:M    
%             X_A_POST(:,ant_idx,user_cnt) = squeeze(X_A_PRI(:,ant_idx,user_cnt)) + ...
%                                            V_A_PRI(ant_idx,user_cnt)/(V_A_PRI(ant_idx,user_cnt)+sigma) * ...
%                                            (A_k')*(Y_k(:,ant_idx)-A_k*squeeze(X_A_PRI(:,ant_idx,user_cnt)));
                                       
            X_A_POST(:,ant_idx,user_cnt) = ((1/sigma)*A_k'*A_k + (1/V_A_PRI(ant_idx,user_cnt))*eye(size(A_k'*A_k)))^(-1) ...
                                        *((1/sigma)*A_k'*Y_k(:,ant_idx) + (1/V_A_PRI(ant_idx,user_cnt))*squeeze(X_A_PRI(:,ant_idx,user_cnt)));
%             
%             V_A_POST(ant_idx,user_cnt) = V_A_PRI(ant_idx,user_cnt)-(Tp/N)*(V_A_PRI(ant_idx,user_cnt))^2/(V_A_PRI(ant_idx,user_cnt)+sigma);
%             if(V_A_POST(ant_idx,user_cnt) <= 0)
%                 V_A_POST(ant_idx,user_cnt) = 1e-20;
%                 info_disp = sprintf('%d-th user Posterior Var of module A is less than 0\n',user_cnt);
%                 disp(info_disp);
%             end
%             V_A_EXT(ant_idx,user_cnt) = min(1e20, max(1e-20, (1/V_A_POST(ant_idx,user_cnt)-1/V_A_PRI(ant_idx,user_cnt))^-1));


            V_temp_1 = ((1/V_A_PRI(ant_idx,user_cnt))/N)*trace(((1/sigma)*A_k'*A_k + (1/V_A_PRI(ant_idx,user_cnt))*eye(size(A_k'*A_k)))^(-1));
            V_temp_2 = (1/V_A_PRI(ant_idx,user_cnt))./V_temp_1;
            V_A_POST(ant_idx,user_cnt) = 1/V_temp_2;
            if(V_A_POST(ant_idx,user_cnt) < 0)
                V_A_POST(ant_idx,user_cnt) = 1e-8;
            end
            V_A_EXT(ant_idx,user_cnt) = min(1e20, max(1e-20, (1/V_A_POST(ant_idx,user_cnt)-1/V_A_PRI(ant_idx,user_cnt))^-1));
            
            X_A_EXT(:,ant_idx,user_cnt) = V_A_EXT(ant_idx,user_cnt) * ...
                                          (squeeze(X_A_POST(:,ant_idx,user_cnt))/V_A_POST(ant_idx,user_cnt) - ...
                                           squeeze(X_A_PRI(:,ant_idx,user_cnt))/V_A_PRI(ant_idx,user_cnt));

            V_B_PRI(ant_idx,user_cnt) = beta.* V_A_EXT(ant_idx,user_cnt) + (1-beta).* V_A_EXT_Last(ant_idx,user_cnt);
            X_B_PRI(:,ant_idx,user_cnt) = beta.* X_A_EXT(:,ant_idx,user_cnt) + (1-beta).* X_A_EXT_Last(:,ant_idx,user_cnt);
            % input message to channel support, alpha = 1 means non-zero
            % coefficient, alpha = -1 means zero-element
            pi_in_alpha(:,ant_idx,user_cnt) = Alpha_input_message(squeeze(X_B_PRI(:,ant_idx,user_cnt)),V_B_PRI(ant_idx,user_cnt),sysPar(user_cnt),ant_idx);
        end
        
    end
    
    V_A_EXT_Last = V_A_EXT;
    X_A_EXT_Last = X_A_EXT;
    
    % Consider binary MRF prior of alpha, and it is a 1st order MRF, which means
    % the MRF is a lattice
    
    % Model the input messages from neighborhoods as Bernoulli distribution
    if(debug_mode)
        parfor user_cnt = 1:K
            supp_matrix = sysPar(user_cnt).supp_mat;
            pi_out_alpha_vec = zeros(N*M,1);
            idx1 = find(reshape(supp_matrix,[],1)==1);
            idx2 = find(reshape(supp_matrix,[],1)==0);
            pi_out_alpha_vec(idx1) = 0.99;
            pi_out_alpha_vec(idx2) = 0.01;
            pi_out_alpha(:,:,user_cnt) = reshape(pi_out_alpha_vec,N,M);
        end
        DebugData(iter_cnt).pi_in_alpha = pi_in_alpha;
    else
        if(Prior_Model == 'HMRF')
            qT = sysParComm.qT;
            [pi_out_alpha, qT_upt, MRF_updates]  = HiddenMRF_Prior_Est(pi_in_alpha, a, b, MRF_Inner_Iteration, qT, EM_learning);
            DebugData(iter_cnt).pi_out_alpha = pi_out_alpha;
            DebugData(iter_cnt).pi_in_alpha = pi_in_alpha;
            DebugData(iter_cnt).MRF_updates = MRF_updates;
            sysParComm.qT = qT_upt;
            a = MRF_updates(2);
            b = MRF_updates(1);
        end
        if(Prior_Model == 'CMRF')
            [pi_out_alpha, S_POST] = MRF_3D_Prior_Est(pi_in_alpha,a,b,c,MRF_Inner_Iteration);
            %% EM update MRF Parameters
            if(EM_learning)
            % self-designed EM method. 
                S_HAT = double(S_POST > 1/2);
                options = optimset('GradObj', 'on', 'Algorithm', 'interior-point', ...
                        'MaxFunEvals', 20,'Display', 'off');
                lb = [0; -1; -1];   % Lower bounds [alpha,beta,gamma]
                ub = [1; 1; 1];    % Upper bounds [alpha,beta,gamma]
                [MRF_updates] = fmincon(@pseudoLF_3D, [a; b; c], [], [], [], ...
                                    [], lb, ub, [], options, S_HAT, N, M, K);
            % EM method from Vila
%                 S_HAT = double(S_POST > 1/2);
%             
%                 PI_pad = 1/2*ones(N+2,M+2,K+2);	% Pad in all dimensions with dummy nodes
% 
%                 PI_pad(2:N+1,2:M+1,2:K+1) = S_POST;  % Padded cube, w/ 1/2 at dummy nodes
% 
%                 ShiftPI_pad = 2*PI_pad - 1;     % 2*pi - 1
% 
%                 NeighborhoodSum = ShiftPI_pad(1:N,2:M+1,2:K+1) + ...
%                                     ShiftPI_pad(3:N+2,2:M+1,2:K+1) + ...
%                                     ShiftPI_pad(2:N+1,1:M,2:K+1) + ...
%                                     ShiftPI_pad(2:N+1,3:M+2,2:K+1) + ...
%                                     ShiftPI_pad(2:N+1,2:M+1,1:K) + ...
%                                     ShiftPI_pad(2:N+1,2:M+1,3:K+2);      % \sum_Neigh(n) (2*pi(q) - 1)
%                 % \sum_Neigh(n) (2*pi(n) - 1) (2*pi(q) - 1)
%                 AvgNeighborhoodSum = (2*S_POST - 1) .* NeighborhoodSum;
% 
%                 options = optimset('GradObj', 'on', 'Hessian', ...
%                         'off', 'MaxFunEvals', 100, 'tolX', 1e-8, 'Display', ...
%                         'notify');
%                 lb_3d = [-1; -3];   % Lower bounds [alpha; beta]
%                 ub_3d = [1; 3];    % Upper bounds [alpha; beta]
% 
%                 [MRF_updates, ~, ~] = fmincon(@meanfieldLF, ...
%                                           [a; b], [], [], [], [], lb_3d, ub_3d, [], ...
%                                           options, S_POST, NeighborhoodSum, AvgNeighborhoodSum);
            else
                MRF_updates = [a;b;c];
            end
            a = MRF_updates(1);
            b = MRF_updates(2);
            c = MRF_updates(3);
            DebugData(iter_cnt).pi_out_alpha = pi_out_alpha;
            DebugData(iter_cnt).pi_in_alpha = pi_in_alpha;
            DebugData(iter_cnt).MRF_updates = MRF_updates;
        end
    end
    
    
    for user_cnt = 1:K
        for ant_idx = 1:M
            column_sparsity = pi_out_alpha(:,ant_idx,user_cnt);
            [sysPar(user_cnt).theta(ant_idx), sysPar(user_cnt).sigVar(ant_idx)] = X_Mean_Var_Est_MRF_B(N,squeeze(X_B_PRI(:,ant_idx,user_cnt)),V_B_PRI(ant_idx,user_cnt), column_sparsity ,sysPar(user_cnt).theta(ant_idx),sysPar(user_cnt).sigVar(ant_idx));
            [X_B_POST_temp, Var_B_temp] = CBG_Est_MRF(X_B_PRI(:,ant_idx,user_cnt),V_B_PRI(ant_idx,user_cnt), column_sparsity ,sysPar(user_cnt).theta(ant_idx),sysPar(user_cnt).sigVar(ant_idx));
            
            X_B_POST(:,ant_idx,user_cnt) = X_B_POST_temp;
            Var_B(:,ant_idx,user_cnt) = Var_B_temp;
            
            V_B_POST(ant_idx,user_cnt) = (mean(squeeze(Var_B(:,ant_idx,user_cnt))));
            
            if(V_B_POST(ant_idx,user_cnt) < 0)
                V_B_POST(ant_idx,user_cnt) = 1e-8;
                info_disp = sprintf('%d-th user Posterior Var of module B is less than 0\n',user_cnt);
                disp(info_disp);
            end
            
            V_B_EXT(ant_idx,user_cnt) = 1/(1/V_B_POST(ant_idx,user_cnt) -1/V_B_PRI(ant_idx,user_cnt));
            
            X_B_EXT(:,ant_idx,user_cnt) = (V_B_EXT(ant_idx,user_cnt)/V_B_POST(ant_idx,user_cnt)) .* squeeze(X_B_POST(:,ant_idx,user_cnt)) - ...
                                          (V_B_EXT(ant_idx,user_cnt)/V_B_PRI(ant_idx,user_cnt)) .* squeeze(X_B_PRI(:,ant_idx,user_cnt));

            V_A_PRI(ant_idx,user_cnt) = beta.* V_B_EXT(ant_idx,user_cnt) + (1-beta).* V_B_EXT_Last(ant_idx,user_cnt);
            X_A_PRI(:,ant_idx,user_cnt) = beta.* X_B_EXT(:,ant_idx,user_cnt) + (1-beta).* X_B_EXT_Last(:,ant_idx,user_cnt);

        end
        V_B_EXT_Last = V_B_EXT;
        X_B_EXT_Last = X_B_EXT;
    end
    
    X_Est = X_B_POST;
    NMSE(iter_cnt) = immse(X_Est,X_org)./immse(X_org,zeros(size(X_org)));
    parfor user_cnt = 1:K
        NMSE_individual(iter_cnt,user_cnt) = immse(squeeze(X_Est(:,:,user_cnt)),squeeze(X_org(:,:,user_cnt)))./immse(squeeze(X_org(:,:,user_cnt)),zeros(size(squeeze(X_org(:,:,user_cnt)))));
    end
end


end