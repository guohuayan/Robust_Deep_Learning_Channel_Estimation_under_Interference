function [pi_out, T_POST, rho10_update, qT_upt] = ModuleB_MU_HMM_Prior(pi_in,qT,lambda,rho10)
    [~,K] = size(pi_in);
    if(K>1)
        parfor k = 1:K
            pi_in_k = pi_in(:,k);
            qT_k = qT(k);
            pi_in_MM_k(:,k) = (pi_in_k.*qT_k + (1-pi_in_k).*(1-qT_k)) ./ ...
                        (pi_in_k.*qT_k + (1-pi_in_k).*(1-qT_k) + (1-pi_in_k));
        end
        pi_in_MM = squeeze(prod(pi_in_MM_k,2)./(prod(pi_in_MM_k,2) + prod((1-pi_in_MM_k),2)));

        [pi_out_MM, T_POST, rho10_update] = ModuleB_HMMPrior(pi_in_MM,lambda,rho10);

        for k = 1:K
            prod_idx = setdiff(1:K,k);
            pi_in_MM_S = pi_in_MM_k(:,k);
            pi_out_MM_k = (pi_out_MM)./(pi_out_MM + (1-pi_out_MM).*prod(((1./squeeze(pi_in_MM_k(:,prod_idx)))-1),2));
            qT_k = qT(k);
            [qT_k_upt]=qT_EM_Update(T_POST, pi_out_MM_k, pi_in_MM_S, qT_k);

            pi_out(:,k) = pi_out_MM_k.* qT_k_upt;
            qT_upt(k) = qT_k_upt;
        end
    else
        pi_in_MM = pi_in;
        [pi_out_MM, T_POST, rho10_update] = ModuleB_HMMPrior(pi_in_MM,lambda,rho10);
        pi_out = pi_out_MM;
        qT_upt = qT;
    end
%     lambda = min(max(0, mean(T_POST)), 1); 
end

