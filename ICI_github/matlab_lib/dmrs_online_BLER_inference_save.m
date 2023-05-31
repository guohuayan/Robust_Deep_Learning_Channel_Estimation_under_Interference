function flag = dmrs_online_BLER_inference_save(speed,SNRdB,He_freq_w,Iot,H_delay_w)
    savename=strcat('dmrs_final_vbion_channel_test_',num2str(speed),'kmh_',num2str(SNRdB),'_iot_10.mat');
    save(savename,'He_freq_w','speed','SNRdB','Iot','H_delay_w');
    flag=1;
end

