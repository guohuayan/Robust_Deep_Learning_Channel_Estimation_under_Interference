function A = regular_smdf(A)
    d=diag(A);
    A=A-d+real(d);
    A=1/2*(A+A');
end

