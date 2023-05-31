function [power_d_o,cov_H_mix_o,He_delay_buff,Cd_int_buff,io] = online_VBI_training_load(speed,SNRdB)
    loadname=strcat('Training_BLER_',num2str(speed),'kmh_snr_',num2str(SNRdB),'pilot_2_iot_10_v3.mat');
    % load(loadname,'power_d_o','cov_H_mix_o','cov_H_mix','power_d','SNRdB','He_delay_buff','Cd_int_buff','io','-v7.3');
    load(loadname,'power_d_o','cov_H_mix_o','He_delay_buff','Cd_int_buff','io');
end

