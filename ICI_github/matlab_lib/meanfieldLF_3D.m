% meanfieldLF function
% *************************************************************
% This is the mean-field approximation of the likelihood
% function
%
% INPUTS:
%   params: A 2-by-1 vector of [beta; alpha] parameter values,
%           at which the EM cost function and gradient are to 
%           be computed
%   PI:     An Nx-by-Ny-by-Nz tensor containing 
%           Pr{s(i,j,k) = 1 | Y}
%   NeighSum:   An Nx-by-Ny-by-Nz tensor whose (i,j,k)th
%               element consists of \sum_q 2*pi(q) - 1, where
%               the sum is over entries, q, that are neighbors
%               of the (i,j,k)th voxel, and pi(q) = 
%               Pr{s(q) = 1 | Y}.
%   AvgNeighSum:    Similar to NeighSum, except that the
%                   (i,j,k)th element is multiplied by
%                   Pr{s(i,j,k) = 1 | Y}.
%   
% OUTPUTS:
%   f:      Cost function at [alpha; beta]
%   g:      Gradient of cost function at [alpha; beta]
%
% Coded by: Justin Ziniel, The Ohio State Univ.
% E-mail: zinielj@ece.osu.edu
% Last change: 04/13/13
% Change summary: 
%       - Created (04/13/13; JAZ)
% Version 0.2

function [f, g] = meanfieldLF_3D(params, PI)
    alpha = params(1);
    beta_x = params(2);
    beta_y = params(3);
    beta_z = params(4);
    
    % *************************************************************
    % Compute the quantities that will be used to compute a
    % mean-field approximation to an expectation-maximization (EM)
    % update of the MRF parameters, alpha and beta.  In what
    % follows, we assume S(n) is drawn from {-1,1}
    [Nx,Ny,Nz]=size(PI);
    PI_pad = 1/2*ones(Nx+2,Ny+2,Nz+2);	% Pad in all dimensions with dummy nodes
%             PI = msgProds(:,:,:,2) ./ sumMsgProds;    % Pr{S(n) = 1 | Y}
%             PI(missingIdx) = 1/2;
    PI_pad(2:Nx+1,2:Ny+1,2:Nz+1) = PI;  % Padded cube, w/ 1/2 at dummy nodes
    PI_pad(1,2:Ny+1,2:Nz+1) = PI(end,:,:);
    PI_pad(end,2:Ny+1,2:Nz+1) = PI(1,:,:);
    PI_pad(2:Nx+1,1,2:Nz+1)= PI(:,end,:);
    PI_pad(2:Nx+1,end,2:Nz+1)= PI(:,1,:);
    PI_pad(2:Nx+1,2:Ny+1,1)= PI(:,:,end);
    PI_pad(2:Nx+1,2:Ny+1,end)= PI(:,:,1);

    % Compute a couple different sums that appear often in the
    % mean-field update expressions
    ShiftPI_pad = 2*PI_pad - 1;     % 2*pi - 1
    NeighSum_x = ShiftPI_pad(1:Nx,2:Ny+1,2:Nz+1) + ...
        ShiftPI_pad(3:Nx+2,2:Ny+1,2:Nz+1) ;
    NeighSum_y = ShiftPI_pad(2:Nx+1,1:Ny,2:Nz+1) + ...
        ShiftPI_pad(2:Nx+1,3:Ny+2,2:Nz+1);
    NeighSum_z = ShiftPI_pad(2:Nx+1,2:Ny+1,1:Nz) + ...
        ShiftPI_pad(2:Nx+1,2:Ny+1,3:Nz+2);      % \sum_Neigh(n) (2*pi(q) - 1)
    % \sum_Neigh(n) (2*pi(n) - 1) (2*pi(q) - 1)
    AvgNeighSum_x = (2*PI - 1) .* NeighSum_x;
    AvgNeighSum_y = (2*PI - 1) .* NeighSum_y;
    AvgNeighSum_z = (2*PI - 1) .* NeighSum_z;

    % Start by computing the objective function value, which is
    % the posterior expectation of the mean field approximation
    % of the 3D MRF prior
    PosSum = -alpha + beta_x*NeighSum_x+beta_y*NeighSum_y+beta_z*NeighSum_z;    % s(n) = 1
    ExpPosSum = exp(PosSum);
    ExpNegSum = exp(-PosSum);
    f = alpha*(2*PI - 1) - beta_x*AvgNeighSum_x- beta_y*AvgNeighSum_y- beta_z*AvgNeighSum_z + log(ExpPosSum + ExpNegSum);
    f = sum(f(:));

    % Next, compute the derivative of the cost function w.r.t.
    % alpha and beta
    g_alpha = (ExpNegSum - ExpPosSum) ./ (ExpNegSum + ExpPosSum);
    g_alpha = (2*PI - 1) + g_alpha;     % deriv of f wrt alpha
    g_beta_x = (NeighSum_x.*ExpPosSum - NeighSum_x.*ExpNegSum) ./ ...
        (ExpNegSum + ExpPosSum);
    g_beta_x = -AvgNeighSum_x + g_beta_x;     % deriv of f wrt beta
    g_beta_y = (NeighSum_y.*ExpPosSum - NeighSum_y.*ExpNegSum) ./ ...
        (ExpNegSum + ExpPosSum);
    g_beta_y = -AvgNeighSum_y + g_beta_y;     % deriv of f wrt beta
    g_beta_z = (NeighSum_z.*ExpPosSum - NeighSum_z.*ExpNegSum) ./ ...
        (ExpNegSum + ExpPosSum);
    g_beta_z = -AvgNeighSum_z + g_beta_z;     % deriv of f wrt beta

    g = [sum(g_alpha(:)); sum(g_beta_x(:)); sum(g_beta_y(:)); sum(g_beta_z(:))];
end