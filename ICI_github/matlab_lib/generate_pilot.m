function [Home_pilot,Intercell_pilot] = generate_pilot(num_pilot,Np,Ci)
    F_pilot=generate_ZC_pilot(Np,1,0);
    Nd=Np*2;
    %%
    Home_pilot=zeros(Nd,num_pilot);
    Home_pilot=generate_rand_QPSK(Home_pilot);
    Home_pilot(1:2:end,:)=repmat(F_pilot,[1,num_pilot]);
    %%
%     Ci=length(Intercell_channel_G);
    %% interference cell
    Intercell_pilot=cell(Ci,1);
    for c0=1:Ci
        %% pilot
        F_pilot_int_base=generate_ZC_pilot(Np,c0+2,0);
        int_pilot=zeros(Nd,num_pilot);
        int_pilot=generate_rand_QPSK(int_pilot);
        int_pilot(1:2:end,:)=repmat(F_pilot_int_base,[1,num_pilot]);
        Intercell_pilot{c0}.F_pilot_int=int_pilot;
    end
end

