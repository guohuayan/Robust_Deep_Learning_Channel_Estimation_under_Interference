function y = generate_rand_QPSK(x)
    X=size(x);
    re=1-2.*randi([0,1],X);
    im=1-2.*randi([0,1],X);
    y=re+1j.*im;
end

