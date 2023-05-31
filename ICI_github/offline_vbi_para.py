import torch
import numpy as np
import scipy.io as io
import mat73
from scipy.io import savemat
from scipy.linalg import dft

from vbi_fun import *


PATH_para = './ICI_github/offline_vbitrain_major_sub60_1e5.mat'
Fm = mat73.loadmat(PATH_para)['Fm']  # load from .mat
Fm=torch.from_numpy(Fm).cfloat()
Fa = mat73.loadmat(PATH_para)['Fa']  # load from .mat
Fa=torch.from_numpy(Fa).cfloat()
Cs = mat73.loadmat(PATH_para)['cov_H_mix']  # load from .mat
Cs=torch.from_numpy(Cs).cfloat()
power_d = mat73.loadmat(PATH_para)['power_d']  # load from .mat
power_d=torch.from_numpy(power_d).cfloat()


PATH_cnn = './ICI_github/major_channel_cdl_off_rddelaypre_sys_para.pt'
torch.save({'Fm':Fm,
            'Fa':Fa,
            'Cs':Cs,'power_d':power_d,
            }, PATH_cnn)

