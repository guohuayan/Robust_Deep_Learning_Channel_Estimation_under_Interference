function [out,csi,nvar] = myEqualizeMMSE_ICI_common_angular(rxSym,Hest,Cn)
%%  rxSym: 2304x64
%%  Hest: 2304x64
%% Cn: 64x64

    % Extract input dimensions
    nRE = size(Hest,1); % Number of resource elements 2304
    R = size(Hest,2);   % Number of receive antennas 64
    P = size(Hest,3);   % Number of transmission planes 1

    % Initialize output based on the input dimensions
    out = zeros(nRE,P,'like',rxSym);
    csi = zeros(nRE,P,class(rxSym));
    %%
    % Permute to R-by-P-NRE to allow quicker access
    Hest = permute(Hest,[2 3 1]); %% 2304 64 1 -> 64   1  2304
    %%
    nvar=real(mean(diag(Cn)));

%     if (R == 1 && P == 1)
%         % For SISO case
%         Hdash = conj(Hest);
%         csi = Hdash.*Hest + nVar;
%         G = Hdash./csi;
%         out = G.*rxSym;
%     else
        % For non-SISO case

        

%         n0eye = Cn;
        invC=inv(Cn);
        for REIdx=1:nRE
            H = Hest(:,:,REIdx);           
            invH = inv(H'*invC*H+eye(P));
            csi(REIdx, :)  = 1./real(diag(invH))*nvar;
            G = invH*H'*invC;
            out(REIdx,:) = G*(rxSym(REIdx,:).');
            %%
%             temp = (H'*H + n0eye)\Itx;
%             csi(REIdx,:) = 1./real(diag(temp));
%             G = temp*H';
%             out(REIdx,:) = G*(rxSym(REIdx,:).');
        end
%     end

end
