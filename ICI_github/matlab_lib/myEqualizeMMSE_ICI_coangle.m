function [Decov_w,csi,nvar] = myEqualizeMMSE_ICI_coangle(H_home,C_int)
%% H_home: 192x64
%% C_int: 64x64x192
%% Cn: 64x64
%% int_power: 1x1

    % Extract input dimensions
    nRE = size(H_home,1); % Number of resource elements 2304
    R = size(H_home,2);   % Number of receive antennas 64
    P = 1;

    %%
    % Permute to R-by-P-NRE to allow quicker access
    H_home=H_home.'; %% 64 192
    %%
    Decov_w=zeros(P,R,nRE);  %% 1,64,192
    csi=zeros(nRE,P);
    for REIdx=1:nRE
        H = H_home(:,REIdx);
        Ci= C_int;
        invC=inv(Ci);
        invH = inv(H'*invC*H+eye(P));
%         invH = inv(H'/(Cn+Ci)*H+eye(P));
        nvar=real(mean(diag(Ci)));
        csi(REIdx, :)  = 1./real(diag(invH))*nvar;
        G = invH*H'*invC;
        Decov_w(:,:,REIdx)=G;
    end
%     end

end
