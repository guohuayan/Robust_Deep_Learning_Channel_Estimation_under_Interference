function [pi_out, T_POST, rho10_update] = my_HMM_Prior(pi_in,lambda,rho10)
%     [~,K] = size(pi_in);
%     K=1;
    pi_in_MM = pi_in;
    [pi_out_MM, T_POST, rho10_update] = Module_HMMPrior(pi_in_MM,lambda,rho10);
    pi_out = pi_out_MM;
end

function [pi_out, T_POST, rho10_update] = Module_HMMPrior(pi_in,lambda,rho10)

    N = length(pi_in);
    Rho_10 = rho10*ones(N,1);
    if lambda==1
        Rho_01=1*ones(N,1);
    else
        Rho_01 = max(min(Rho_10 .* (lambda ./ (1 - lambda)),1),0);
    end
    gamma_FWD = NaN*ones(N,1);
    gamma_BWD = NaN*ones(N,1);
%     gamma_FWD(1,:) = lambda;
%     gamma_BWD(N,:) = 0.5;
    gamma_FWD(1,:) = eps;
    gamma_BWD(N,:) = eps;
    % First the forward pass
    for n = 2 : 1 : N
        FW_num = (Rho_01(n,:) .* (1 - pi_in(n-1,:)) .* (1 - gamma_FWD(n-1,:)) + (1 - Rho_10(n,:)) .* pi_in(n-1,:) .* gamma_FWD(n-1,:));
        FW_den =  ((1 - pi_in(n-1,:)) .* (1 - gamma_FWD(n-1,:)) + pi_in(n-1,:) .* gamma_FWD(n-1,:));
        gamma_FWD(n,:) = FW_num./FW_den;
    end
    % Now the backward pass
    for n = N-1 : -1 : 1
        BW_num = (Rho_10(n,:).* (1-gamma_BWD(n+1,:)) .* (1-pi_in(n+1,:)) + (1 - Rho_10(n,:)) .* pi_in(n+1,:) .* gamma_BWD(n+1,:));
        BW_den = BW_num + (Rho_01(n,:).* pi_in(n+1,:) .* gamma_BWD(n+1,:) + (1-Rho_01(n,:)).*(1-pi_in(n+1,:)).* (1-gamma_BWD(n+1,:)));
        gamma_BWD(n,:) = BW_num./BW_den;
    end
    
    pi_out = (gamma_FWD .* gamma_BWD) ./ ((1 - gamma_FWD) .* (1 - gamma_BWD) + gamma_FWD .* gamma_BWD);
    %% EM Estimation Needs here
    % First compute posterior means
    % E = E(P(t))*R(P(s))
    Mean_T = (pi_in .* gamma_FWD .* gamma_BWD) ./ ...
             ((1 - pi_in) .* (1 - gamma_FWD) .* (1 - gamma_BWD) + ...
             pi_in .* gamma_FWD .* gamma_BWD); 

    P_0_0_Factoral = (1 - Rho_01(1:N-1,:)) .* ...
                     ((1 - gamma_FWD(1:N-1,:)) .* (1 - pi_in(1:N-1,:))) .* ...
                     ((1 - gamma_BWD(2:N,:)) .* (1 - pi_in(2:N,:)));

    P_0_1_Factoral = Rho_01(1:N-1,:) .* ...
                    ((1 - gamma_FWD(1:N-1,:)) .* (1 - pi_in(1:N-1,:))) .* ...
                    ((gamma_BWD(2:N,:)) .* (pi_in(2:N,:)));


    P_1_0_Factoral = Rho_10(1:N-1,:) .* ...
                    ((gamma_FWD(1:N-1,:)) .* (pi_in(1:N-1,:))) .* ...
                    ((1 - gamma_BWD(2:N,:)) .* (1 - pi_in(2:N,:)));

    P_1_1_Factoral = (1 - Rho_10(1:N-1,:)) .* ...
                     ((gamma_FWD(1:N-1,:)) .* (pi_in(1:N-1,:))) .* ...
                     ((gamma_BWD(2:N,:)) .* (pi_in(2:N,:)));

    E_corr = (P_1_1_Factoral)./((P_0_0_Factoral + P_0_1_Factoral + P_1_0_Factoral + P_1_1_Factoral));

    Rho_10_update = sum(sum(Mean_T(1:N-1) - E_corr))./sum(sum(Mean_T(1:N-1)));

    rho10_update = max(min(Rho_10_update, 1), 0);

    T_POST = Mean_T;

end


