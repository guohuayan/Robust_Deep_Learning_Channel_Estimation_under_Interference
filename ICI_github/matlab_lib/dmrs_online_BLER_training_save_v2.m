function flag = dmrs_online_BLER_training_save_v2(speed,power_d_o,cov_H_mix_o,cov_H_mix,power_d,SNRdB,He_delay_buff,Cd_int_buff,io)
    savename=strcat('dmrs_BLER_training_',num2str(speed),'kmh_snr_',num2str(SNRdB),'_iot_10_v2.mat');
    save(savename,'power_d_o','cov_H_mix_o','cov_H_mix','power_d','SNRdB','He_delay_buff','Cd_int_buff','io','-v7.3');
    flag=1;
end

