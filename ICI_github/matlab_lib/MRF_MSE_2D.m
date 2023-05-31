function [theta_update,phi_update,a_upt,b_upt,X_B_EXT,X_B_POST,V_B_EXT] = MRF_MSE_2D(X_B_PRI,V_B_PRI,phi,theta,a,b,InnerIter,EM_Update,BG_Mode)
% phi,%Variance of the signal
% theta, %Mean of signal
% a, MRF parameter alpha, refer to local correlation parameter.
% b, MRF parameter beta, refer to cross correlation parameter.
% InnerIter, MRF inner iteration
% EM_Update, set to 1 enable MRF parameter update.
% BG_Mode, disable MRF prior.
[N,M] = size(X_B_PRI);

X_B_PRI_Vec = reshape(X_B_PRI,[],1);

tmpVar = V_B_PRI + phi;

pi_in=(tmpVar./V_B_PRI).*exp((abs(X_B_PRI_Vec - theta).^2 ./tmpVar)-(abs(X_B_PRI_Vec).^2 ./V_B_PRI));

pi_alpha_out = 1./ (1 + pi_in);

if(BG_Mode == 0)
    [pi_out_alpha,S_POST] = MRF_Prior_Est(reshape(pi_alpha_out,size(X_B_PRI)),a,b,InnerIter);

    %% EM update MRF Parameters
    if(EM_Update)
        S_HAT = double(S_POST > 1/2);

        options = optimset('GradObj', 'on', 'Algorithm', 'interior-point', ...
                'MaxFunEvals', 20,'Display', 'off');
        lb = [-1; 0];   % Lower bounds [beta; alpha]
        ub = [1; 1];    % Upper bounds [beta; alpha]
        [MRF_updates] = fmincon(@pseudoLF, [b; a], [], [], [], ...
                [], lb, ub, [], options, S_HAT, N, M);
        a_upt = MRF_updates(2);
        b_upt = MRF_updates(1);
    else
        a_upt = a;% larger a encourage sparser matrix
        b_upt = b;% larger b encourage larger non-zero clusters. 
    end
else
    pi_out_alpha = reshape(pi_alpha_out,size(X_B_PRI));
    a_upt = a;% larger a encourage sparser matrix
    b_upt = b;% larger b encourage larger non-zero clusters. 
end
[theta_update, phi_update] = X_Mean_Var_Est_MRF_B(N*M,X_B_PRI_Vec,V_B_PRI,reshape(pi_out_alpha,[],1),theta,phi);
[X_B_POST_Vec, Var_B] = CBG_Est_MRF(X_B_PRI_Vec,V_B_PRI, reshape(pi_out_alpha,[],1) ,theta_update,phi_update);

X_B_POST = reshape(X_B_POST_Vec,size(X_B_PRI));

V_B_POST = mean(Var_B);
V_B_EXT = 1/(1/V_B_POST -1/V_B_PRI);
X_B_EXT_Vec = (V_B_EXT/V_B_POST) * X_B_POST_Vec- (V_B_EXT/V_B_PRI) * X_B_PRI_Vec;

X_B_EXT = reshape(X_B_EXT_Vec,size(X_B_PRI));
end

