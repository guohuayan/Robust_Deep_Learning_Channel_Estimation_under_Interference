close all
clear all

load('major_test_CDLA_sub60_5e2_1.mat','H_home_batch_tmp','H_int_batch_tmp');
[Nfft,M,len]=size(H_home_batch_tmp);
Nslot=len*6;
H_home_batch=zeros(Nfft,M,Nslot);
H_int_batch=zeros(Nfft,M,2,Nslot);
idx=1:len;
H_home_batch(:,:,idx)=H_home_batch_tmp;
H_int_batch(:,:,:,idx)=H_int_batch_tmp;
%%
load('major_test_CDLA_sub60_5e2_2.mat','H_home_batch_tmp','H_int_batch_tmp');
idx=idx+len;
H_home_batch(:,:,idx)=H_home_batch_tmp;
H_int_batch(:,:,:,idx)=H_int_batch_tmp;
load('major_test_CDLA_sub60_5e2_3.mat','H_home_batch_tmp','H_int_batch_tmp');
idx=idx+len;
H_home_batch(:,:,idx)=H_home_batch_tmp;
H_int_batch(:,:,:,idx)=H_int_batch_tmp;
load('major_test_CDLA_sub60_5e2_4.mat','H_home_batch_tmp','H_int_batch_tmp');
idx=idx+len;
H_home_batch(:,:,idx)=H_home_batch_tmp;
H_int_batch(:,:,:,idx)=H_int_batch_tmp;
load('major_test_CDLA_sub60_5e2_5.mat','H_home_batch_tmp','H_int_batch_tmp');
idx=idx+len;
H_home_batch(:,:,idx)=H_home_batch_tmp;
H_int_batch(:,:,:,idx)=H_int_batch_tmp;
load('major_test_CDLA_sub60_5e2_6.mat','H_home_batch_tmp','H_int_batch_tmp');
idx=idx+len;
H_home_batch(:,:,idx)=H_home_batch_tmp;
H_int_batch(:,:,:,idx)=H_int_batch_tmp;
%%
save('major_test_CDLA_sub60_3e3.mat','H_home_batch','H_int_batch');






