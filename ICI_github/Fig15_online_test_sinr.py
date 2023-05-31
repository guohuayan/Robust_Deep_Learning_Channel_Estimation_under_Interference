# The prior model in this version will be more simple
import os
from pickle import TRUE
import mat73
import scipy.io as io

# from vbi_fun import *
from online_fine_fun import *
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

class CustomOnlineDataset(Dataset):
    def __init__(self, H_input, mse,  H_label):
        self.H_input=H_input
        self.mse=mse
        self.H_delay=H_label
        self.len=H_label.shape[0]
        self.shape=H_label.shape
        # self.transform = transform
        # self.target_transform = target_transform

    def __len__(self):
        return self.len

    def __getitem__(self, idx):
        input = self.H_input[idx]
        mse = self.mse[idx]
        label = self.H_delay[idx]
        return input, mse, label

    # to obtain data sample number, Nc and Nt
    def get_dataset_shape(self):
        return self.shape

def test_loop_module_B(model,dataloader,  loss_fn, Fm,Fa, device = torch.device('cpu')):
    Fm=Fm.to(device)
    Fa=Fa.to(device)
    # mse_est=mse_est.to(device)
    num_batches = len(dataloader)
    test_loss = 0
    model.eval()
    # pr_w = 0
    model.set_to_inference()

    with torch.no_grad():
        for i, data in enumerate(dataloader, 0):
            # get the inputs; data is a list of [inputs, labels]
            inputs = data[0].to(device, dtype=torch.cfloat)
            err = data[1].to(device, dtype=torch.float)
            labels = data[2].to(device, dtype=torch.cfloat)
            # outputs = model.reconstruct_img(inputs)
            # outputs, z, b =model(inputs)

            # err=torch.square( torch.abs(inputs - labels)).mean(dim=(1,2))
            # err=mse_est.expand(inputs.shape[0])
            # err=mse_est

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

class OnlineDataset(Dataset):
    def __init__(self, HLS_delay_w,Ct_int_w, batch_len):
        self.H_delay=HLS_delay_w
        self.Int_Cov=Ct_int_w
        # self.Cu_inv=Cu_inv_w

        self.len=self.H_delay.shape[0]
        self.batch_len=batch_len

        self.cursor=0

    def __len__(self):
        return self.len

    def __getitem__(self):
        idx=torch.arange(self.batch_len)+self.cursor
        H_delay = self.H_delay[idx]
        Int_Cov = self.Int_Cov[idx]
        # Cu_inv = self.Cu_inv[idx]
        self.cursor=self.cursor+1
        return H_delay,  Int_Cov

    # to obtain data sample number, Nc and Nt
    def reset_cursor(self):
        self.cursor=0
        return None

def pre_process_mat(x,shift):
    x=torch.from_numpy(x)
    if shift==1:
        x=x.permute(2,0,1); # 1000 128 32
    x=x.cfloat()          
    return x

def check_test_point(Cs_test,power_d_test,Ht_delay_test,Cn,C_int_test,conFa,H_test,
                vae,test_criterion,Fm,Fa):
    Ht_vbi, mse_est= LMMSE_vbi(Ht_delay_test,Cs_test,power_d_test,Cn,C_int_test)

    Ht_est = torch.matmul(Fm,Ht_vbi)
    vbi_err=NMSE(Ht_est, H_test)

    Ht_vbi=torch.matmul(Ht_vbi,conFa)
    test_dataset =CustomOnlineDataset(H_input=Ht_vbi, mse=mse_est, H_label=H_test)
    testloader = DataLoader(dataset=test_dataset, batch_size=128, 
                                    shuffle=False, num_workers=1, pin_memory=True,drop_last=True)
    test_loss=test_loop_module_B(vae, testloader, test_criterion,Fm,Fa, device=device)          
    return test_loss, vbi_err

if __name__ == '__main__':
    device = torch.device('cuda:3' if torch.cuda.is_available() else 'cpu')
    # device = 'cpu'
    print(device)
    print('Using device:', device)

    # PATH_train='./no_GB/save_result/dmrs_channel_pre_cnn_1e5.pt'
    # PATH_train='./major/save_data/dmrs_channel_3gpp_pre_cnn_1e5.pt'
    # PATH_train= './major/save_data/major_channel_cdl_off_rddelaypre_cnn_1e5.pt'
    PATH_train= './ICI_github/major_channel_cdl_off_rddelaypre_sys_para.pt'
    

    Fm=torch.load(PATH_train)['Fm']
    conFm=torch.conj(Fm.permute(1,0))
    Fa=torch.load(PATH_train)['Fa']
    conFa=torch.conj(Fa.permute(1,0))

    PATH_tpara = './ICI_github/offline_vbitrain_major_sub60_1e5.mat'
    Cs_off = mat73.loadmat(PATH_tpara)['cov_H_mix']  # load from .mat
    Cs_off=torch.from_numpy(Cs_off).cfloat()
    power_d_off = mat73.loadmat(PATH_tpara)['power_d']  # load from .mat
    power_d_off =torch.from_numpy(power_d_off).cfloat()

    PATH_tpara = './ICI_github/major_online_vbitrain_CDLA_sub60_test_1e5.mat'
    Cs_true = mat73.loadmat(PATH_tpara)['cov_H_mix']  # load from .mat
    Cs_true=torch.from_numpy(Cs_true).cfloat()
    power_d_true = mat73.loadmat(PATH_tpara)['power_d']  # load from .mat
    power_d_true =torch.from_numpy(power_d_true).cfloat()

    n_ant=Cs_off.shape[0]
    n_delay=power_d_off.shape[0]

    
    Model_PATH = "./ICI_github/save_result/major_cdlmodel_rddelay_noici_cbn_Filtervbi_para_n20_3_cdlA_v3.pt"
    # Above is best
    # Model_PATH = "./major/save_result/major_cdlmodel_rddelay_noici_cbn_Filtervbi_para_n20_5_cdlA_v4.pt"
    #  not good in low snr
    # Model_PATH = "./major/save_result/major_cdlmodel_rddelay_noici_cbn_Filtervbi_para_n20_3_cdlA_v3_tn10_0.pt"
    # slightly worse than the best  v3 n20_3
    # Model_PATH = "./major/save_result/robust_noici_CNN_model_rddelay_hyper_para_n20_5_cdlA_testn10_0_t1.pt"
    #  not good in low snr
    # Model_PATH = "./major/save_result/major_cdlmodel_rddelay_noici_cbn_Filtervbi_para_n20_3_cdlA_v3_meanall.pt"
    # almost the same as the best
    # Model_PATH = "./major/save_result/major_cdlmodel_rddelay_noici_cbn_Filtervbi_para_n20_5_cdlA_v4.pt"
     # better in high snr but worse in low snr
    
    vae = VAE(device=device)
    vae.load_state_dict(torch.load(Model_PATH))
    # vae.to(device)
    for name, param in vae.named_parameters():
        param.requires_grad = False
    vae.to(device)

    snr_w=torch.linspace(-10,0,5)
    mse_vbioff=torch.zeros_like(snr_w)
    mse_nnoff=torch.zeros_like(snr_w)
    mse_vbigenie=torch.zeros_like(snr_w)
    mse_nngenie=torch.zeros_like(snr_w)


    for i, test_inr in enumerate(snr_w):
        if i==0:
            PATH_data='./ICI_github/online_test_CDLA_60k_sinr_n10_Inr_3.mat'
        elif i==1:
            PATH_data='./ICI_github/online_test_CDLA_60k_sinr_n7.5_Inr_3.mat'
        elif i==2:
            PATH_data='./ICI_github/online_test_CDLA_60k_sinr_n5_Inr_3.mat'
        elif i==3:
            PATH_data='./ICI_github/online_test_CDLA_60k_sinr_n2.5_Inr_3.mat'
        elif i==4:
            PATH_data='./ICI_github/online_test_CDLA_60k_sinr_n0_Inr_3.mat'

        HLS_delay_w = mat73.loadmat(PATH_data)['HLS_delay_w']  # Np,M,ite
        HLS_delay_w = pre_process_mat(HLS_delay_w,1)
        
        Ct_int_w= mat73.loadmat(PATH_data)['Ct_int_w']  # Np,M,ite
        Ct_int_w = pre_process_mat(Ct_int_w,1)

        # Cu_inv_w= mat73.loadmat(PATH_data)['Cu_inv_w']  # Np,M,ite
        # Cu_inv_w = pre_process_mat(Cu_inv_w,0)

        noise_p = mat73.loadmat(PATH_data)['noise_p']  # Np,M,ite
        noise_p = pre_process_mat(noise_p,0)

        Ht_freq_w= mat73.loadmat(PATH_data)['Ht_freq_w']  # Np,M,ite
        Ht_freq_w = pre_process_mat(Ht_freq_w,1)

        # online_dataset =OnlineDataset(HLS_delay_w=HLS_delay_w[1000:3000,:,:],
        #                     Ct_int_w=Ct_int_w[1000:3000,:,:],batch_len=20)
        # online_dataset.reset_cursor()
    
        Cn=noise_p*torch.eye(n_ant,dtype=torch.cfloat)

        H_test = Ht_freq_w[0:1000,:,:]
        Ht_delay_test=HLS_delay_w[0:1000,:,:]
        C_int_test=Ct_int_w[0:1000,:,:]

        Cs_test=Cs_off
        power_d_test=power_d_off    
        test_loss, vbi_err=check_test_point(Cs_test,power_d_test,Ht_delay_test,Cn,C_int_test,conFa,H_test,
                    vae,test_criterion,Fm,Fa)
        print(vbi_err) 
        print(test_loss) 
        mse_vbioff[i]=vbi_err
        mse_nnoff[i]=test_loss 

        Cs_test=Cs_true 
        power_d_test=power_d_true     
        test_loss, vbi_err=check_test_point(Cs_test,power_d_test,Ht_delay_test,Cn,C_int_test,conFa,H_test,
                    vae,test_criterion,Fm,Fa)
        print(vbi_err) 
        print(test_loss) 
        mse_vbigenie[i]=vbi_err
        mse_nngenie[i]=test_loss   

    mse_vbioff=mse_vbioff.numpy()
    mse_nnoff=mse_nnoff.numpy()
    mse_vbigenie=mse_vbigenie.numpy()
    mse_nngenie=mse_nngenie.numpy()

    snr_w=snr_w.numpy()
    print("Done test!")
    print(mse_nngenie)
    mdic = {"mse_vbioff": mse_vbioff,"mse_nnoff":mse_nnoff,
            "mse_vbigenie": mse_vbigenie,"mse_nngenie":mse_nngenie,
                "snr_w":snr_w}   

    savemat("./ICI_github/CNN_major_snrvary_online_ici_test_cdlA_varysinr_inr_3.mat", mdic) 
    # savemat("./major/save_result/CNN_major_snrvary_online_ici_test_cdlA_varysinr_inr_3_test.mat", mdic)     

    
    


     