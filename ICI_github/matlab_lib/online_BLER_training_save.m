function flag = online_BLER_training_save(speed,SNRdB,power_d_o,cov_H_mix_o,cov_H_mix,power_d,He_delay_buff,Cd_int_buff,io)
    savename=strcat('Training_BLER_',num2str(speed),'kmh_snr_',num2str(SNRdB),'pilot_2_iot_10_v3.mat');
    save(savename,'power_d_o','cov_H_mix_o','cov_H_mix','power_d','SNRdB','He_delay_buff','Cd_int_buff','io','-v7.3');
    flag=1;
end

