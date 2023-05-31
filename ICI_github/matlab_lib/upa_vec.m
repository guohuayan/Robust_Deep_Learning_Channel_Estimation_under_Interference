function alpha = upa_vec(N1,N2,phi,xi)
    phi = wrapAzimuthAngles(phi); % wrap azimuth angles to [-180,180]
    xi = wrapZenithAngles(xi); % wrap zenith angles to [0,360] and map [180,360] to [180,0]
    phi=phi/180*pi;
    xi=xi/180*pi;
    alpha1=1:N1;
    alpha1=alpha1-1;
    alpha1=1j.*pi.*sin(xi).*sin(phi).*alpha1;
    alpha1=exp(alpha1);
    alpha1=alpha1.';
    alpha2=-(1:N2);
    alpha2=alpha2+1;
    alpha2=1j.*pi.*cos(xi).*alpha2;
    alpha2=exp(alpha2);
    alpha2=alpha2.';
    alpha=kron(alpha1,alpha2);
end

