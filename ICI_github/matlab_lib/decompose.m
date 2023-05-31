function U = decompose(A)
    [V,D] = eig(A); %% A = V*D*V'.
    s=diag(D);
    s=real(s);
    s=sqrt(s);
    s=max(s,0);
    s=diag(s);
    U=V*s;
end

