# The prior model in this version will be more simple
import os
from pickle import TRUE

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

def LMMSE_ici_online(x,Cs,alpha,Cn,Int_Cov,device = torch.device('cpu')):
    x_shape = x.shape  # 1000 15 64
    batch_len=x_shape[0]   
    n_delay = x_shape[1]
    n_ant =  x_shape[2]
    Cs_w=torch.empty(n_delay,n_ant,n_ant,dtype=torch.cfloat,device =device)
    for i in torch.arange(n_delay):
        Cs_w[i,:,:]=alpha[i]*Cs  # 400 64 64

    H_delay=x
    x_est=torch.zeros_like(H_delay)
    C_err=torch.zeros(n_delay,n_ant,n_ant)
    for i in torch.arange(batch_len):
        Ci=Int_Cov[i,:,:].reshape(1,n_ant,n_ant)
        Ci=Ci.expand(n_delay,n_ant,n_ant)
        ts=Cs_w+Cn.expand(n_delay,n_ant,n_ant)+Ci
        tmp_inv=torch.linalg.inv(ts)
        filter=torch.matmul(Cs_w,tmp_inv)
        h_delay = H_delay[i,:,:].reshape(n_delay,n_ant,1)
        x_est[i,:,:]=torch.matmul(filter,h_delay).reshape(n_delay,n_ant)

        ti=torch.linalg.inv(Cs_w)+torch.linalg.inv(Cn.expand(n_delay,n_ant,n_ant)+Ci)
        C_err=C_err+torch.linalg.inv(ti)
    C_err=C_err/batch_len

    return x_est, C_err

def LMMSE_ici_online_speed(x,Cs,alpha,Cn,Int_Cov,device = torch.device('cpu')):
    x_shape = x.shape  # 1000 15 64
    batch_len=x_shape[0]   
    n_delay = x_shape[1]
    n_ant =  x_shape[2]
    Cs_inv=torch.linalg.inv(Cs)
    Cs_w_inv=torch.empty(n_delay,n_ant,n_ant,dtype=torch.cfloat,device =device)
    for i in torch.arange(n_delay):
        Cs_w_inv[i,:,:]=Cs_inv/alpha[i]  # 400 64 64

    H_delay=x
    x_est=torch.zeros_like(H_delay)
    C_err=torch.zeros(n_delay,n_ant,n_ant)
    for i in torch.arange(batch_len):
        Ci_inv=torch.linalg.inv(Int_Cov[i,:,:]+Cn).reshape(1,n_ant,n_ant)
        Ci_inv=Ci_inv.expand(n_delay,n_ant,n_ant)
        ti=torch.linalg.inv(Cs_w_inv+Ci_inv)
        filter=torch.matmul(ti,Ci_inv)

        # Ci=Int_Cov[i,:,:].reshape(1,n_ant,n_ant)
        # Ci=Ci.expand(n_delay,n_ant,n_ant)
        # ts=Cs_w+Cn.expand(n_delay,n_ant,n_ant)+Ci
        # tmp_inv=torch.linalg.inv(ts)
        # filter=torch.matmul(Cs_w,tmp_inv)
        h_delay = H_delay[i,:,:].reshape(n_delay,n_ant,1)
        x_est[i,:,:]=torch.matmul(filter,h_delay).reshape(n_delay,n_ant)

        # ti=torch.linalg.inv(Cs_w)+torch.linalg.inv(Cn.expand(n_delay,n_ant,n_ant)+Ci)
        # C_err=C_err+torch.linalg.inv(ti)
        C_err=C_err+ti
    C_err=C_err/batch_len

    return x_est, C_err

def LMMSE_nature_grad(power_d,Cs,Covx,rho,device = torch.device('cpu')):
    n_delay=Covx.shape[0]
    n_ant=Covx.shape[1]

    Covi=torch.linalg.inv(Cs)
    Covi=Covi.expand(n_delay,n_ant,n_ant)
    power=torch.matmul(Covx,Covi)
    npower_d=torch.diagonal(power,dim1=1, dim2=2).mean(-1)

    ext_power=npower_d.reshape(n_delay,1,1)
    ext_power=ext_power.expand(n_delay,n_ant,n_ant)

    # Cxe=torch.div(Covx,ext_power.real).mean(dim=0)
    Cxe=((torch.tensor(1.0,dtype=torch.float,device =device)/ext_power.real)*Covx).mean(dim=0)

    nCs=Cxe/torch.trace(Cxe)*n_ant

    dpower_d=npower_d-power_d
    dCs=torch.linalg.inv(nCs)-torch.linalg.inv(Cs)

    power_d=power_d+rho*dpower_d
    Cs=torch.linalg.inv(Cs)+rho*dCs
    Cs=torch.linalg.inv(Cs)

    return power_d, Cs

def LMMSE_est_Cs(x, C_err,device = torch.device('cpu')):
    x_shape = x.shape
    batch_len=x_shape[0]   
    n_delay = x_shape[1]
    n_ant =  x_shape[2]

    C_h=torch.zeros(n_delay,n_ant,n_ant)
    for i in torch.arange(batch_len):
        h=x[i,:,:].reshape(n_delay,n_ant,1)
        conh=torch.conj(h.permute(0,2,1))
        C_ht=torch.matmul(h,conh)
        C_h=C_h+C_ht
    C_h=C_h/batch_len

    Covx=C_h+C_err

    return Covx

# This function is wrong when alpha is very small
def LMMSE_vbi_speed(x,Cs,alpha,Cn,Int_Cov,device = torch.device('cpu')):
    x_shape = x.shape  # 1000 15 64   
    n_delay = x_shape[1]
    n_ant =  x_shape[2]
    Cs_inv=torch.linalg.inv(Cs)
    Cs_w_inv=torch.empty(n_delay,n_ant,n_ant,dtype=torch.cfloat,device =device)
    for i in torch.arange(n_delay):
        Cs_w_inv[i,:,:]=(torch.tensor(1.0,dtype=torch.float,device =device)/alpha[i].real)*Cs_inv  # 400 64 64

    H_delay=x
    x_est=torch.zeros_like(H_delay)
    batch_len=x.shape[0]
    est_mse_w=torch.zeros(batch_len,1,dtype=torch.float,device =device)
    for i in torch.arange(batch_len):
        Ci=Int_Cov[i,:,:]
        Cu=Ci+Cn
        Cu_inv=torch.linalg.inv(Cu)
        est_tmp=torch.zeros(n_delay,dtype=torch.cfloat,device =device)
        for j in torch.arange(n_delay):
            Cs_i=Cs_w_inv[j,:,:]
            C_ebs=torch.linalg.inv(Cu_inv+Cs_i)
            filter=torch.matmul(C_ebs,Cu_inv)
            h_delay = H_delay[i,j,:].reshape(n_ant,1)
            x_est[i,j,:]=torch.matmul(filter,h_delay).reshape(n_ant)
            est_tmp[j]=torch.trace(C_ebs)/n_ant
        est_mse_w[i]=est_tmp.mean().real

    return x_est,est_mse_w

def LMMSE_vbi(x,Cs,alpha,Cn,Int_Cov,device = torch.device('cpu')):
    x_shape = x.shape  # 1000 15 64
    batch_len=x_shape[0]   
    n_delay = x_shape[1]
    n_ant =  x_shape[2]
    Cs_w=torch.empty(n_delay,n_ant,n_ant,dtype=torch.cfloat,device =device)
    for i in torch.arange(n_delay):
        Cs_w[i,:,:]=alpha[i]*Cs  # 400 64 64

    H_delay=x
    x_est=torch.zeros_like(H_delay)
    batch_len=x.shape[0]
    est_mse_w=torch.zeros(batch_len,1,dtype=torch.float,device =device)
    alpha_thr=torch.tensor(1e-8,dtype=torch.float,device =device)
    Cs_inv=torch.linalg.inv(Cs)
    for i in torch.arange(batch_len):
        Ci=Int_Cov[i,:,:].reshape(1,n_ant,n_ant)
        Ci=Ci.expand(n_delay,n_ant,n_ant)
        ts=Cs_w+Cn.expand(n_delay,n_ant,n_ant)+Ci
        tmp_inv=torch.linalg.inv(ts)
        filter=torch.matmul(Cs_w,tmp_inv)
        h_delay = H_delay[i,:,:].reshape(n_delay,n_ant,1)
        x_est[i,:,:]=torch.matmul(filter,h_delay).reshape(n_delay,n_ant)

        est_tmp=torch.zeros(n_delay,dtype=torch.cfloat,device =device)
        Cu_inv=torch.linalg.inv(Cn+Int_Cov[i,:,:])
        for j in torch.arange(n_delay):
            if alpha[j].real>alpha_thr:
                Cs_i=(torch.tensor(1.0,dtype=torch.float,device =device)/alpha[j].real)*Cs_inv
                C_ebs=torch.linalg.inv(Cu_inv+Cs_i)
                est_tmp[j]=torch.trace(C_ebs)/n_ant

        est_mse_w[i]=est_tmp.mean().real

        # ti=torch.linalg.inv(Cs_w)+torch.linalg.inv(Cn.expand(n_delay,n_ant,n_ant)+Ci)
        # C_err=torch.linalg.inv(ti)
        # est_mse_w[i]=torch.diagonal(C_err,0,1,2).mean().real
    return x_est, est_mse_w

def LMMSE(x,Cs,alpha,Cn,device = torch.device('cpu')):
    x_shape = x.shape  # 1000 15 64   
    n_delay = x_shape[1]
    n_ant =  x_shape[2]
    Filter=torch.empty(n_delay,n_ant,n_ant,dtype=torch.cfloat,device =device)
    for i in torch.arange(n_delay):
        tmp=alpha[i]*Cs  # 400 64 64
        ts=tmp+Cn
        tmp1=torch.linalg.inv(ts)
        tmp=torch.matmul(tmp,tmp1)
        Filter[i,:,:]=tmp
    
    # x = x.reshape(x_shape[0],x_shape[1],x_shape[2],1)
    # Filter=Filter.reshape(1,n_delay,64,64)
    # x_est = torch.matmul(Filter,x)
    # x_est = x_est.reshape(x_shape)

    H_delay=x
    x_est=torch.zeros_like(H_delay)
    batch_len=x.shape[0]
    for i in torch.arange(batch_len):
        h_delay=torch.squeeze(H_delay[i,:,:]).reshape(n_delay,n_ant,1)
        h_est=torch.matmul(Filter,h_delay).reshape(n_delay,n_ant)
        x_est[i,:,:]=h_est

    return x_est


def trans_to_torch(x, iscomplex=True):
    x = torch.from_numpy(x)
    if iscomplex:
        x = x.cfloat()  
    else:
        x = x.float()       
    return x

def NMSE(output, label):
    output=torch.view_as_real(output)
    label=torch.view_as_real(label)
    loss1 = torch.square(output - label).sum(dim=(1,2,3))
    loss2 = torch.square( label).sum(dim=(1,2,3))
    loss = torch.div(loss1,loss2).mean()           
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
        
        
    