function xest = amp(yn,var,A,niter)
% support: the length of nonzero element

[M, N] = size(A);
delta = M/N;

lambda=sqrt(log(N)*var);
% lambda=sqrt(var);

yp = yn;
z = yp;
tau = 1/M*(norm(z,2)^2);

%% Approximate Message Passing for basis selection
% Initializing
xest = zeros(N,1);

for t=1:niter
    
    % depending on if there is noise in the measuerements (amp0) the denoiser
    % parameter is calculated differenetly
    deniser_parameter = lambda+tau;
    
    eta_prime = 1/N * nnz(wthresh(A'*z + xest,'s', deniser_parameter));
    
    xest = wthresh(A'*z + xest,'s', deniser_parameter);
%     sum(xest)
    
    z = yp - A*xest + 1/delta * z * eta_prime;
    
    tau=eta_prime/delta*(lambda+tau);
    
%     z = yp - A*xest + 1/delta * z * eta_prime;
    

%% the original way to get the denoiser parameter
%     eta_prime = 1/N * nnz(wthresh(A'*z(:,t-1) + x(:,t),'s', deniser_parameter));   
%     if(amp0 == true)
%         tau2(t) = tau2(t-1)/delta * eta_prime;
%     else
%         tau2(t) = (1+tau2(t-1))/delta * eta_prime;
%     end
%%    
%     tau2(t) = 1/M*(norm(z(:,t),2)^2);

   
end

% for iter=1:niter
%    sigma = norm(z,2)/sqrt(M);
%    xest = wthresh(xest + A'*z,'s', tau*sigma);
% %    xest = eta(xest + phi'*z, tau*sigma);
%    tmp = nnz(xest );
%    z = yp - A*xest + (tmp/M)*z;
% end

    
end
