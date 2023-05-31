import torch
import numpy as np
import scipy.io as io
import mat73
from scipy.io import savemat
from scipy.linalg import dft

from vbi_fun import *

# PATH_train='./major/save_data/dmrs_channel_3gpp_pre_cnn_1e5.pt'
PATH_train= './major/save_data/major_channel_cdl_off_rddelaypre_cnn_1e5.pt'
# H_home=torch.load(PATH_train)['H_home'] # 1000 128 32

Fm=torch.load(PATH_train)['Fm']
conFm=torch.conj(Fm.permute(1,0))
Fa=torch.load(PATH_train)['Fa']
# Fa=torch.kron(torch.eye(2),Fa.contiguous())
conFa=torch.conj(Fa.permute(1,0))
Cs=torch.load(PATH_train)['Cs']
power_d=torch.load(PATH_train)['power_d']


# snr_w=torch.linspace(-20,10,10000)
snr_w=torch.linspace(-20,5,10000)
# snr_w=torch.linspace(-15,5,10000)
Filter_w=LMMSE_filter(snr_w,Cs,power_d)

# PATH_para = './no_GB/save_data/offline_vbitrain_1e5.mat'
# PATH_para = './major/matlabcode/offline_vbitrain_3gpp_1e5.mat'
PATH_para = './major/matlabcode/offline_vbitrain_major_sub60_1e5.mat'
R_w = mat73.loadmat(PATH_para)['R_w']  # load from .mat
R_w=torch.from_numpy(R_w).cfloat()
R_w=R_w.permute(2,0,1)
n_ant=R_w.shape[-1]
n_delay=R_w.shape[0]

R_inv=torch.zeros_like(R_w)
for i in torch.arange(n_delay):
    R_tmp=R_w[i,:,:]
    R_inv[i,:,:]=torch.linalg.inv(R_tmp)

mse_est_w=torch.zeros_like(snr_w)
for i, snr in enumerate(snr_w):
    inv_noise_var=torch.float_power(torch.tensor(10), snr/torch.tensor(10))
    inv_noise_var=inv_noise_var.cfloat()
    Cn_inv=inv_noise_var*torch.eye(n_ant,dtype=torch.cfloat)
    mse_est=torch.empty(n_delay,1,dtype=torch.float)
    for j in torch.arange(n_delay):
        ts=R_inv[j,:,:]+Cn_inv
        tmp1=torch.linalg.inv(ts)
        tmp=torch.trace(tmp1)/n_ant
        mse_est[j]=tmp.real
    mse_est_w[i]=mse_est.mean()
mse_inv=1/mse_est_w
weights=mse_inv/mse_inv.mean()
print(weights.min())
print(weights.max())
# print(weights)

# PATH_cnn = './major/save_result/LMSE_filter_rddlay_snr_n15_5.pt'
# PATH_cnn = './major/save_result/LMSE_filter_rddlay_snr_n20_10.pt'
# PATH_cnn = './major/save_result/LMSE_filter_rddlay_snr_n20_3.pt'
PATH_cnn = './major/save_result/LMSE_filter_rddlay_snr_n20_5.pt'
torch.save({'snr_w':snr_w,
            'Filter_w':Filter_w,
            'Fm':Fm,
            'Fa':Fa,
            'weights':weights,'R_inv':R_inv,
            'Cs':Cs,'power_d':power_d,
            }, PATH_cnn)

