# The prior model in this version will be more simple
import os
from pickle import TRUE
import ssl

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

def LMMSE_ici(x,Cs,alpha,Cn,Int_Cov,device = torch.device('cpu')):
    x_shape = x.shape  # 1000 15 64   
    n_delay = x_shape[-2]
    n_ant =  x_shape[-1]
    Cs_w=torch.empty(n_delay,n_ant,n_ant,dtype=torch.cfloat,device =device)
    for i in torch.arange(n_delay):
        Cs_w[i,:,:]=alpha[i]*Cs  # 400 64 64

    H_delay=x
    x_est=torch.zeros_like(H_delay)
    if x.dim()==2:
        Ci=Int_Cov.reshape(1,n_ant,n_ant)
        Ci=Ci.expand(n_delay,n_ant,n_ant)
        ts=Cs_w+Cn.expand(n_delay,n_ant,n_ant)+Ci
        tmp_inv=torch.linalg.inv(ts)
        filter=torch.matmul(Cs_w,tmp_inv)
        h_delay = H_delay.reshape(n_delay,n_ant,1)
        x_est=torch.matmul(filter,h_delay).reshape(n_delay,n_ant)
    else:
        batch_len=x.shape[0]
        for i in torch.arange(batch_len):
            Ci=Int_Cov[i,:,:].reshape(1,n_ant,n_ant)
            Ci=Ci.expand(n_delay,n_ant,n_ant)
            ts=Cs_w+Cn.expand(n_delay,n_ant,n_ant)+Ci
            tmp_inv=torch.linalg.inv(ts)
            filter=torch.matmul(Cs_w,tmp_inv)
            h_delay = H_delay[i,:,:].reshape(n_delay,n_ant,1)
            x_est[i,:,:]=torch.matmul(filter,h_delay).reshape(n_delay,n_ant)

    return x_est

def LMMSE(x,Cs,alpha,Cn,device = torch.device('cpu')):
    x_shape = x.shape  # 1000 15 64   
    n_delay = x_shape[-2]
    n_ant =  x_shape[-1]
    Filter=torch.empty(n_delay,n_ant,n_ant,dtype=torch.cfloat,device =device)
    for i in torch.arange(n_delay):
        tmp=alpha[i]*Cs  # 400 64 64
        ts=tmp+Cn
        tmp1=torch.linalg.inv(ts)
        tmp=torch.matmul(tmp,tmp1)
        Filter[i,:,:]=tmp

    H_delay=x
    x_est=torch.zeros_like(H_delay)
    if x.dim()==2:
        h_delay=torch.squeeze(H_delay).reshape(n_delay,n_ant,1)
        h_est=torch.matmul(Filter,h_delay).reshape(n_delay,n_ant)
        x_est=h_est
    else:
        batch_len=x.shape[0]
        for i in torch.arange(batch_len):
            h_delay=torch.squeeze(H_delay[i,:,:]).reshape(n_delay,n_ant,1)
            h_est=torch.matmul(Filter,h_delay).reshape(n_delay,n_ant)
            x_est[i,:,:]=h_est

    return x_est

def LMMSE_fast(x,Filter_w,idx_w,device = torch.device('cpu')):
    x_shape = x.shape  # 1000 15 64   
    n_delay = x_shape[-2]
    n_ant =  x_shape[-1]

    H_delay=x
    x_est=torch.zeros_like(H_delay)
    for i, idx in enumerate(idx_w):
        h_delay=torch.squeeze((H_delay[i,:,:])).reshape(n_delay,n_ant,1)
        Filter=torch.squeeze(Filter_w[:,:,:,idx])
        h_est=torch.matmul(Filter,h_delay).reshape(n_delay,n_ant)
        x_est[i,:,:]=h_est

    return x_est

def LMMSE_filter(snr_w,Cs,alpha,device = torch.device('cpu')):   
    n_delay = alpha.shape[0]
    n_ant =  Cs.shape[0]
    n_snr=snr_w.shape[0]

    Filter_w=torch.empty(n_delay,n_ant,n_ant,n_snr,dtype=torch.cfloat,device =device)
    for s0, snr in enumerate(snr_w):
        noise_var=torch.float_power(torch.tensor(10,device=device), -snr/torch.tensor(10,device=device))
        noise_var=noise_var.cfloat()
        Cn=noise_var*torch.eye(n_ant,dtype=torch.cfloat,device=device)
        Filter=torch.empty(n_delay,n_ant,n_ant,dtype=torch.cfloat,device =device)
        for i in torch.arange(n_delay):
            tmp=alpha[i]*Cs  # 400 64 64
            ts=tmp+Cn
            tmp1=torch.linalg.inv(ts)
            tmp=torch.matmul(tmp,tmp1)
            Filter[i,:,:]=tmp
        Filter_w[:,:,:,s0]=Filter

    return Filter_w

def NMSE(output, label):
    output=torch.view_as_real(output)
    label=torch.view_as_real(label)
    loss1 = torch.square(output - label).sum(dim=(1,2,3))
    loss2 = torch.square( label).sum(dim=(1,2,3))
    loss = torch.div(loss1,loss2).mean()           
    return loss

def NMSE_element(output, label):
    output=torch.view_as_real(output)
    label=torch.view_as_real(label)
    loss1 = torch.square(output - label).sum(dim=(1,2,3))
    loss2 = torch.square( label).sum(dim=(1,2,3))
    loss = torch.div(loss1,loss2)          
    return loss

def dft_column(N, x):
    """
    A column of DFT matrix
    :param N: size of DFT matrix
    :param x:  value of the column
    :return:    a column corresponding to x
    """
    k = np.arange(N)
    f = np.exp(- 1j * 2 * np.pi * x * k)

    return f


def DFTMatrix(N):
    """
    Return the DFT matrix F in Dai's paper
    :param N:   Size of DFT matrix -- N x N
    :return:    F -- the DFT matrix F in Dai's paper
    """
    # construct the DFT matrix
    F = np.zeros((N, N), dtype=complex)
    for n in range(N):
        x = n / N
        F[:, n] = dft_column(N, x)

    return F
        
        
    