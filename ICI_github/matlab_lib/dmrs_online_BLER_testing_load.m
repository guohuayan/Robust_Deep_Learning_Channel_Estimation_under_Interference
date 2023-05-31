function channel_sample = dmrs_online_BLER_testing_load(speed,SNRdB)
    loadname=strcat('dmrs_final_channel_test_',num2str(speed),'kmh_',num2str(SNRdB),'_iot_10.mat');
    load(loadname,'channel_sample');
end

