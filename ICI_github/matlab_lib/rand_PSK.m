function y = rand_PSK(size)
    phi=randn(size)*2*pi;
    y=exp(1j.*phi);
%     y=sqrt(1/2).*(randn(size)+1j.*randn(size));
end