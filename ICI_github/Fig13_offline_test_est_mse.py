# The input mse is estimated instead of the true mse
import os
from pickle import TRUE
import mat73
import scipy.io as io

from vbi_fun import *
from major_CNN_noise_noici_cbn_hyper_fun import *

import numpy as np
import torch
import torch.nn as nn
import torchvision.transforms as transforms
from torch.distributions import constraints
from torch.utils.data import Dataset,DataLoader
import torch.nn.functional as F
import torch.optim as optim


import torchvision
import matplotlib.pyplot as plt
from scipy.io import savemat



def test_loop(model,dataloader,  loss_fn, Fm,Fa, mse_est,device = torch.device('cpu')):
    Fm=Fm.to(device)
    Fa=Fa.to(device)
    mse_est=mse_est.to(device)
    num_batches = len(dataloader)
    test_loss = 0
    model.eval()
    # pr_w = 0
    model.set_to_inference()

    with torch.no_grad():
        for i, data in enumerate(dataloader, 0):
            # get the inputs; data is a list of [inputs, labels]
            inputs, labels = data[0].to(device, dtype=torch.cfloat), data[1].to(device, dtype=torch.cfloat)
            # outputs = model.reconstruct_img(inputs)
            # outputs, z, b =model(inputs)

            # err=torch.square( torch.abs(inputs - labels)).mean(dim=(1,2))
            err=mse_est.expand(inputs.shape[0])

            outputs, z, b =model(inputs,err)
            
            outputs=post_NN_module(outputs,Fm,Fa)
            loss = loss_fn(outputs, labels)
            test_loss += loss.item()
            # pr_w = pr_w+ torch.mean(b)

    # test_loss = test_loss/size/M/K
    test_loss = test_loss/num_batches
    # pr_w = pr_w/num_batches
    
    return test_loss

def test_criterion(output, target):
    # loss = torch.square(output - target).sum()
    output=torch.view_as_real(output)
    target=torch.view_as_real(target)
    loss1 = torch.square(output - target).sum(dim=(1,2,3))
    loss2 = torch.square( target).sum(dim=(1,2,3))
    loss = torch.div(loss1,loss2).mean()           
    return loss

if __name__ == '__main__':
    device = torch.device('cuda:2' if torch.cuda.is_available() else 'cpu')
    # device = 'cpu'
    print(device)
    print('Using device:', device)
    
    # PATH_train='./no_GB/save_result/dmrs_channel_pre_cnn_1e5.pt'
    # PATH_train='./major/save_data/dmrs_channel_3gpp_pre_cnn_1e5.pt'
    PATH_train= './ICI_github/major_channel_cdl_off_rddelaypre_sys_para.pt'

    Fm=torch.load(PATH_train)['Fm']
    conFm=torch.conj(Fm.permute(1,0))
    Fa=torch.load(PATH_train)['Fa']
    conFa=torch.conj(Fa.permute(1,0))
    # Cs=torch.load(PATH_train)['Cs']
    # power_d=torch.load(PATH_train)['power_d']

    PATH_tpara = './ICI_github/offline_vbitrain_major_sub60_1e5.mat' 
    Cs = mat73.loadmat(PATH_tpara)['cov_H_mix']  # load from .mat
    Cs=torch.from_numpy(Cs).cfloat()
    power_d = mat73.loadmat(PATH_tpara)['power_d']  # load from .mat
    power_d =torch.from_numpy(power_d).cfloat()


    # H_test=torch.load(PATH_train)['H_test']
    PATH_t ='./ICI_github/major_off_test_mixmodel_sub60_1e3.mat'
    H_test=mat73.loadmat(PATH_t)['H_home_batch']
    H_test=torch.from_numpy(H_test)
    H_test=H_test.permute(2,0,1); # 1000 96 64
    H_test=H_test.cfloat()
    # H_test=H_test_w[0:2000,:,:]

    Ht_delay=torch.matmul(conFm,H_test) # 1000 36 32
    Ht_delay=torch.matmul(Ht_delay,conFa) # 1000 36 32


    
    Model_PATH = "./ICI_github/save_result/major_cdlmodel_rddelay_noici_cbn_Filtervbi_para_n20_3_cdlA_v3.pt"


    # vae = VAE(Fm=Fm, device=device)
    vae = VAE(device=device)
    vae.load_state_dict(torch.load(Model_PATH))
    vae.to(device)

    snr_w=torch.linspace(-10,0,5)
    noise_var_w=-snr_w
    mse_test=noise_var_w
    batch_size_test = 128

    for i, test_snr in enumerate(snr_w):
        Ht_vbi=pre_module(H_test,test_snr,Fm,Cs,power_d)
        Ht_est = torch.matmul(Fm,Ht_vbi)
        err=NMSE(Ht_est, H_test)
        print(err)
        Ht_vbi=torch.matmul(Ht_vbi,conFa)

        test_dataset =CustomImageDataset(H_input=Ht_vbi, H_label=H_test)
        # test_dataset =CustomImageDataset(H_input=Ht_vbi, H_label=Ht_delay)
        testloader = DataLoader(dataset=test_dataset, batch_size=batch_size_test, 
                                    shuffle=False, num_workers=1, pin_memory=True,drop_last=True)

        inv_noise_var=torch.float_power(torch.tensor(10), test_snr/torch.tensor(10))
        inv_noise_var=inv_noise_var.cfloat()
        x_shape = Ht_vbi.shape  # 1000 15 64   
        n_delay = x_shape[-2]
        n_ant =  x_shape[-1]
        Cn_inv=inv_noise_var*torch.eye(n_ant,dtype=torch.cfloat)
        mse_est=torch.empty(n_delay,1,dtype=torch.float)
        for j in torch.arange(n_delay):
            tmp=power_d[j]*Cs  # 400 64 64
            ts=torch.linalg.inv(tmp)+Cn_inv
            tmp1=torch.linalg.inv(ts)
            tmp=torch.trace(tmp1)/n_ant
            mse_est[j]=tmp.real
        mse_est=mse_est.mean()

        test_loss=test_loop(vae, testloader, test_criterion,Fm,Fa, mse_est, device=device)
        print(test_loss)
        mse_test[i]=test_loss


    mse_test=mse_test.numpy()
    snr_w=snr_w.numpy()
    print("Done test!")
    mdic = {"test_loss": mse_test,"snr_w":snr_w}   
    savemat("./ICI_github/major_ici_free_60k_cnn_offtest.mat", mdic)

     