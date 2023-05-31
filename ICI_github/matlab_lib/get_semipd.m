function A = get_semipd(A)
    [V,D] = eig(A); %% A = V*D*V'.
    d=diag(D);
    d=real(d);
    d=max(d,0);
    D=diag(d);
    A = V*D*V';
end

