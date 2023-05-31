import torch
import numpy as np
import scipy.io as io
import mat73
from scipy.io import savemat
from scipy.linalg import dft

from vbi_fun import *

PATH = './ICI_github/dmrs_channel_home_mixmodel_major_sub60_1e5.mat'

H_home = mat73.loadmat(PATH)['H_home_batch']  # Np,M,ite
H_home=torch.from_numpy(H_home)
H_home=H_home.permute(2,0,1); # 1000 128 32
H_home_w=H_home.cfloat()
H_home=H_home_w

PATH_para = './ICI_github/offline_vbitrain_major_sub60_1e5.mat'
Fm = mat73.loadmat(PATH_para)['Fm']  # load from .mat
Fm=torch.from_numpy(Fm).cfloat()
Fa = mat73.loadmat(PATH_para)['Fa']  # load from .mat
Fa=torch.from_numpy(Fa).cfloat()
Cs = mat73.loadmat(PATH_para)['cov_H_mix']  # load from .mat
Cs=torch.from_numpy(Cs).cfloat()
power_d = mat73.loadmat(PATH_para)['power_d']  # load from .mat
power_d=torch.from_numpy(power_d).cfloat()


# PATH_cnn = './major/save_data/major_channel_cdl_off_pre_cnn_1e5.pt'
PATH_cnn = './ICI_github/major_channel_cdl_off_rddelaypre_cnn_1e5.pt'
torch.save({'H_home':H_home,
            'Fm':Fm,
            'Fa':Fa,
            'Cs':Cs,'power_d':power_d,
            }, PATH_cnn)

print(H_home.shape) # 1000 96 64
