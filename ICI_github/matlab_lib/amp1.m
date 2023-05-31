function xest = amp1(yn,sparsity,A,niter)

N = size(A,2);
M = size(A,1);
k=min(sparsity,M);
tau = sqrt(2*log10(M/k)); % needs to be tuned in case of unknown prior sparsity level

%% Approximate Message Passing for basis selection
% Initializing
xest = zeros(N,1);
yp = yn;
z = yp;
eta = @(x,beta) (x./abs(x)).*(abs(x)-beta).*(abs(x)-beta > 0); % denoising function
for iter=1:niter
   sigma = norm(z,2)/sqrt(M);
   xest = eta(xest + A'*z, tau*sigma);
   tmp = sum(abs(xest) > 0);
   z = yp - A*xest + (tmp/M)*z;
   err=yp - A*xest;
    sum(abs(err(:)))
    sum(abs(xest(:)))
end

%[abs(xest) abs(x)]
end