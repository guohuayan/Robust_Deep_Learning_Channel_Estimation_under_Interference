function G = radiation_power_pattern(ZoA,AoA)
    %% vertical
    Av_db=-min(30,12*((ZoA-90)/65)^2);
    %% Horizontal
    Ah_db=-min(30,12*((AoA)/65)^2);
    %%
    A=-min(30,-(Av_db+Ah_db))+8;
    G=10^(A/10);
end