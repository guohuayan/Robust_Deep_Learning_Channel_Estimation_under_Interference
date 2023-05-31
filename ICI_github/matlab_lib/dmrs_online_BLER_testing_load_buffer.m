function [power_d_o,cov_H_mix_o,He_delay_buff,Cd_int_buff,io] = dmrs_online_BLER_testing_load_buffer(speed,SNRdB)
    loadname=strcat('dmrs_BLER_training_',num2str(speed),'kmh_snr_',num2str(SNRdB),'_iot_10.mat');
%              strcat('dmrs_BLER_training_',num2str(speed),'kmh_snr_',num2str(SNRdB),'_iot_10.mat');
    load(loadname,'power_d_o','cov_H_mix_o','He_delay_buff','Cd_int_buff','io');
end

