function [d_2D,h_UT,h_BS] = get_location_config_3gpp(scenario)
    if isequal(scenario,'UMi')
        h_BS=10;
        h_UT=1;
        min_d2D=10;
        ISD=200;
    elseif isequal(scenario,'UMa')
        h_BS=25;
        h_UT=1;
        min_d2D=35;
        ISD=500;
    elseif isequal(scenario,'RMa')
        h_BS=35;
        h_UT=1.5;
        min_d2D=35;
        ISD=1732;
    elseif isequal(scenario,'InH')
        h_BS=3;
        h_UT=1;
        min_d2D=1;
        ISD=20;
    end
    d_2D=(rand(1,1)*0.5+0.5)*ISD/2;
    d_2D=max(d_2D,min_d2D);
end

