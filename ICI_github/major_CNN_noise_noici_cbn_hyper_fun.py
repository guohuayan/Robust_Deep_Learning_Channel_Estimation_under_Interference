# The prior model in this version will be more simple
# redesign the BN
import os
from pickle import TRUE

from vbi_fun import *

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

class CustomImageDataset(Dataset):
    def __init__(self, H_input,  H_label):
        self.H_input=H_input
        # self.H_obs=H_obs
        self.H_delay=H_label
        self.len=H_label.shape[0]
        self.shape=H_label.shape
        # self.transform = transform
        # self.target_transform = target_transform

    def __len__(self):
        return self.len

    def __getitem__(self, idx):
        input = self.H_input[idx]
        # obs = self.H_obs[idx]
        label = self.H_delay[idx]
        return input,  label

    # to obtain data sample number, Nc and Nt
    def get_dataset_shape(self):
        return self.shape


class Decoder(nn.Module):
    def __init__(self, device = torch.device('cpu')):
        super().__init__()
        self.device = device
        # self.Fm=Fm
        # self.Fm=nn.Parameter(self.Fm)
        # self.Fa=Fa
        # self.Fa=nn.Parameter(self.Fa)

        self.conv1 = nn.Sequential(
            nn.Conv2d(2, 16, (3, 5), padding='same'),  # out = (batch, 16, Nc, Nt)
            nn.BatchNorm2d(16),  # out = (batch, 16, Nc, Nt)
            nn.PReLU(16)  # out = (batch, 16, Nc, Nt)
        )

        self.Res_block1=Res_block(device=self.device)

        self.conv4 = nn.Sequential(
            nn.Conv2d(16, 4, (5,5), padding='same'),  # out = (batch, 16, Nc, Nt)
            nn.BatchNorm2d(4),  # out = (batch, 16, Nc, Nt)
            nn.PReLU(4)  # out = (batch, 16, Nc, Nt)
        )

        self.conv5_mean = nn.Sequential(
            nn.Conv2d(4, 2, (5, 5), padding='same'),  # out = (batch, 2, Nc, Nt)
        )

    def forward(self, x,e_out):

        x = torch.view_as_real(x)  # 1000 5 40 2
        x = x.permute(0,3,1,2).contiguous() # 1000 2 15 32
        out = self.conv1(x)

        out = self.Res_block1(out,e_out)

        z = self.conv4(out)
        z = self.conv5_mean(z)
        z = z.permute(0,2,3,1).contiguous() # 1000 15 32 2
        z = torch.view_as_complex(z) # 1000 5 32
        
        # z = torch.matmul(z,self.Fa) # 1000 5 2 32
        # z = z.reshape(z.shape[0],z.shape[1],-1)

        # z = torch.matmul(self.Fm,z) # 1000 96 64
        return z

class prior_z(nn.Module):
    def __init__(self, is_train, device = torch.device('cpu')):
        super().__init__()

        self.is_train = is_train
        self.device = device

        self.conv1 = nn.Sequential(
            nn.Conv2d(2, 4, (3, 3), padding='same'),  # out = (batch, 16, Nc, Nt)
            nn.BatchNorm2d(4),  # out = (batch, 16, Nc, Nt)
            # nn.ReLU(4)  # out = (batch, 16, Nc, Nt)
        )

        self.active1=nn.ReLU()

        self.conv2 = nn.Sequential(
            nn.Conv2d(4, 4, (3, 3), padding='same'),  # out = (batch, 16, Nc, Nt)
            nn.BatchNorm2d(4),  # out = (batch, 16, Nc, Nt)
            # nn.ReLU(4)  # out = (batch, 16, Nc, Nt)
        )

        self.active2=nn.ReLU()

        self.conv3_mean = nn.Sequential(
            nn.Conv2d(4, 1, (3, 3), padding='same'),  # out = (batch, 16, Nc, Nt)
            # nn.Sigmoid(), 
        )

        self.Filter_bank1 = Filter_bank(device=self.device)
        self.Filter_bank2 = Filter_bank(device=self.device)

    def forward(self, x,e_out):
        x = torch.abs(x) # 1000 2 15 32 
        out = self.conv1(x)
        out = self.Filter_bank1(out,e_out)
        out = self.active1(out)


        out = self.conv2(out)
        out = self.Filter_bank2(out,e_out)
        out = self.active2(out)


        out = self.conv3_mean(out) # 1000 1 15 32
        out = torch.squeeze(out) # 1000 15 32

        if self.is_train:
            dM = torch.distributions.gumbel.Gumbel(torch.zeros(out.shape,device=self.device), 
                    torch.ones(out.shape,device=self.device)*torch.tensor(1,device=self.device))
            out = out +dM.sample()
            b=torch.sigmoid(out)
        else:
            b=torch.sigmoid(out)
            b=torch.round(b) 
    
        return b

class Filter_bank(nn.Module):
    def __init__(self, device = torch.device('cpu')):
        super().__init__()
        self.nlinear=nn.Linear(4, 2)

    def forward(self, x, e_out):
        para=self.nlinear(e_out)
        batch_len=para.shape[0]
        gamma=para[:,0].reshape(batch_len,1,1,1)
        beta =para[:,1].reshape(batch_len,1,1,1)
        out=gamma*x+beta
    
        return out

class Res_block(nn.Module):
    def __init__(self, device = torch.device('cpu')):
        super().__init__()
        self.device=device

        self.conv1 = nn.Sequential(
            nn.Conv2d(16, 16, (5,  5), padding='same'),                    # out = (batch, 256, 8-Nc/4, 8-Nt/4)
            nn.BatchNorm2d(16),                                          # out = (batch, 2, Nc, Nt)
        )

        self.active1=nn.ReLU()

        self.conv2 = nn.Sequential(
            nn.Conv2d(16, 16, (5,  5), padding='same'),                    # out = (batch, 256, 8-Nc/4, 8-Nt/4)
            nn.BatchNorm2d(16),                                          # out = (batch, 2, Nc, Nt)
        )

        self.active2=nn.ReLU()

        self.Filter_bank1 = Filter_bank(device=self.device)
        self.Filter_bank2 = Filter_bank(device=self.device)

    def forward(self, x, e_out):

        out1 = self.conv1(x)
        out1 = self.Filter_bank1(out1,e_out)
        out1 = self.active1(out1)

        out2 = self.conv2(out1)
        out2 = self.Filter_bank2(out2,e_out)

        out = out2+x
        out = self.active2(out)
    
        return out

class Encoder(nn.Module):
    def __init__(self, device = torch.device('cpu')):
        super().__init__()
        # setup the three linear transformations used
        self.device=device
        # self.to(device)

        self.conv1 = nn.Sequential(
            nn.Conv2d(2, 16, (9, 7), padding='same'),                    # out = (batch, 256, 8-Nc/4, 8-Nt/4)
            nn.BatchNorm2d(16),                                          # out = (batch, 2, Nc, Nt)
            nn.PReLU(16),                                                # out = (batch, 2, Nc, Nt)
        )

        self.Res_block1=Res_block(device=self.device)

        self.Res_block2=Res_block(device=self.device)

        self.conv3_plus = nn.Sequential(
            nn.Conv2d(16, 16, (3,  3), padding='same'),                    # out = (batch, 256, 8-Nc/4, 8-Nt/4)
            nn.BatchNorm2d(16),                                          # out = (batch, 2, Nc, Nt)
            nn.PReLU(16),                                                # out = (batch, 2, Nc, Nt)
        )

        self.conv4_mean = nn.Sequential(
            nn.Conv2d(16, 2, (3, 3), padding='same'),                     # out = (batch, 2, Nc, Nt)
            # nn.BatchNorm2d(2,affine=False),                                          # out = (batch, 2, Nc, Nt)
        )
        self.conv4_b = nn.Sequential(
            nn.Conv2d(16, 2, (3, 3), padding='same'),                    # out = (batch, 256, Nc, Nt)
            # nn.Sigmoid(),                                               # assume b will not exceed 1
        )       

    def forward(self, x, e_out):

        x = torch.view_as_real(x) # 1000 15 32 -> 1000 15 32 2
        x = x.permute(0,3,1,2).contiguous() # 1000 2 15 32
        out = self.conv1(x)

        out = self.Res_block1(out,e_out)
        out = self.Res_block2(out,e_out)
        out = self.conv3_plus(out)

        z = self.conv4_mean(out)

        # b = z # 1000 2 15 32

        z = z.permute(0,2,3,1).contiguous() # 1000 15 32 2
        z_est = torch.view_as_complex(z) # 1000 15 32

        b = self.conv4_b(out) # 1000 2 15 32
        
        return z_est, b

class VAE(nn.Module):
    # by default our latent space is 50-dimensional
    # and we use 400 hidden units
    def __init__(self,  is_train=True,  device = torch.device('cpu')):
        super().__init__()
        # create the encoder and decoder networks
        self.device=device
        # self.snr = snr
        self.is_train=is_train

        self.encoder = Encoder(device=self.device)
        self.decoder = Decoder(device=self.device)
        self.prior_z = prior_z(self.is_train,device=self.device)

        # self.soft=nn.Softmax(dim=0)

        self.nlinear = nn.Sequential(
            nn.Linear(1, 4),
            nn.ReLU(),
            nn.Linear(4, 4),
            nn.ReLU()
        )      

        self.to(device)

    def forward(self, x, err):
        err=err.reshape(err.shape[0],1)
        e_out=self.nlinear(err)

        x_est = self.encoder(x, e_out)
        z, b = x_est # z: 1000 15 2 32,  b: 1000 4 15 32
        
        H_nn = self.decoder(z,e_out) # 1000 96 64

        b = self.prior_z(b,e_out) # 1000 15 32

        return H_nn, z, b 
    
    def test_speed(self, x, err):
        err=err.reshape(err.shape[0],1)
        e_out=self.nlinear(err)

        x_est = self.encoder(x, e_out)
        z, b = x_est # z: 1000 15 2 32,  b: 1000 4 15 32
        
        H_nn = self.decoder(z,e_out) # 1000 96 64

        # b = self.prior_z(b,e_out) # 1000 15 32

        return H_nn

    def set_to_train(self):
        self.is_train = True

    def set_to_inference(self):
        self.is_train = False

def post_NN_module(z,Fm,Fa):
    z = torch.matmul(z,Fa) # 1000 5 32
    # z = z.reshape(z.shape[0],z.shape[1],-1)

    z = torch.matmul(Fm,z) # 1000 96 64

    return z


def pre_module(H_home,snr,Fm,Cs,power_d,device = torch.device('cpu')):
    noise_var=torch.float_power(torch.tensor(10,device=device), -snr/torch.tensor(10,device=device))
    noise_var=noise_var.cfloat()
    noise=torch.randn(H_home.shape,dtype=torch.cfloat,device=device)*torch.sqrt(noise_var)
    H_LS=H_home+noise
    # H_LS=H_LS.cfloat()

    # n_delay=Fm.size(dim=1)
    conFm=torch.conj(Fm.permute(1,0))
    
    n_ant =  H_home.shape[-1]

    H_delay=torch.matmul(conFm,H_LS) # 1000 15 64
    # delay_shape=H_delay.shape
    # scale=torch.tensor(96/120,device=device)
    Cn=noise_var*torch.eye(n_ant,dtype=torch.cfloat,device=device)
    # Cn=Cn.cfloat()

    H_vbi= LMMSE(H_delay,Cs,power_d,Cn,device)
    # H_est = torch.matmul(Fm,H_vbi)
    # err=NMSE(H_est, H_freq)
    # print(err)
    
    return H_vbi

def Add_noise_hdelay(H_home,snr,Fm,device = torch.device('cpu')):
    noise_var=torch.float_power(torch.tensor(10,device=device), -snr/torch.tensor(10,device=device))
    noise_var=noise_var.cfloat()
    noise=torch.randn(H_home.shape,dtype=torch.cfloat,device=device)*torch.sqrt(noise_var)
    H_LS=H_home+noise
    # H_LS=H_LS.cfloat()

    # n_delay=Fm.size(dim=1)
    conFm=torch.conj(Fm.permute(1,0))
    
    n_ant =  H_home.shape[-1]

    H_delay=torch.matmul(conFm,H_LS) # 1000 15 64

    return H_delay

def pre_module_online(H_home,snr,Fm,Cs,power_d,device = torch.device('cpu')):
    noise_var=torch.float_power(torch.tensor(10,device=device), -snr/torch.tensor(10,device=device))
    noise_var=noise_var.cfloat()
    noise=torch.randn(H_home.shape,dtype=torch.cfloat,device=device)*torch.sqrt(noise_var)
    H_LS=H_home+noise
    # H_LS=H_LS.cfloat()

    # n_delay=Fm.size(dim=1)
    conFm=torch.conj(Fm.permute(1,0))
    
    n_ant =  H_home.shape[-1]

    H_delay=torch.matmul(conFm,H_LS) # 1000 15 64
    # delay_shape=H_delay.shape
    # scale=torch.tensor(96/120,device=device)
    Cn=noise_var*torch.eye(n_ant,dtype=torch.cfloat,device=device)
    # Cn=Cn.cfloat()

    H_vbi= LMMSE(H_delay,Cs,power_d,Cn,device)
    # H_est = torch.matmul(Fm,H_vbi)
    # err=NMSE(H_est, H_freq)
    # print(err)
    
    return H_vbi,H_LS

def est_vbi_para(H_delay,device = torch.device('cpu')):
    h_batch=H_delay # 1000 15 64
    L=h_batch.shape[1]
    M=h_batch.shape[2]
    ite=h_batch.shape[0]
    h_batch_1=h_batch.view(-1,L,M,1)
    h_batch_2=torch.conj(h_batch_1.permute(0,1,3,2))
    R_w=torch.matmul(h_batch_1,h_batch_2).sum(dim=0) # 15 64 64

    Eta=torch.eye(M,device=device)
    alpha=M*ite
    for x in range(100):
        beta_w=torch.zeros(L,device=device,dtype=torch.float)
        for l0 in range(L):
            beta_w[l0]=torch.trace((R_w[l0,:,:].real)*Eta)
        Eta_inv=torch.zeros(M,M,device=device,dtype=torch.cfloat)
        for l0 in range(L):
            Eta_inv=Eta_inv+alpha/beta_w[l0]*R_w[l0,:,:]
        Eta_inv=Eta_inv/L/ite
        ebs=(torch.trace(Eta_inv)).real/M
        Eta_inv=Eta_inv/ebs
        beta_w=beta_w*ebs
        Eta_new=torch.linalg.inv(Eta_inv)
        Eta=Eta_new

    Cs=Eta_inv
    power_d=beta_w/alpha

    return Cs,power_d

def est_vbi_para_eye(H_delay,device = torch.device('cpu')):
    # h_batch=H_delay # 1000 15 64
    # L=h_batch.shape[1]
    M=H_delay.shape[2]
    # ite=h_batch.shape[0]

    obs_sp=torch.view_as_real(H_delay)
    mse1 = torch.square( obs_sp).sum(dim=(3)) # 1000 15 64
    power_d = mse1.mean(dim=(0,2)).cfloat() # 15


    Cs=torch.eye(M,device=device).cfloat()
    # power_d=beta_w/alpha

    return Cs,power_d