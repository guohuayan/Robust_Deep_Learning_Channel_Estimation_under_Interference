# use the true cx for test
import os
from pickle import TRUE

from vbi_fun import *
# from robust_noici_CNN_noise_hyper_fun_v2 import *
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

import mat73



def train(net,trainloader, criterion,optimizer, Fm,Fa,Cs,power_d,snr_w,weights,device = torch.device('cpu')):
    # initialize loss accumulator
    running_loss = 0.0
    Fm=Fm.to(device)
    Fa=Fa.to(device)
    conFm=torch.conj(Fm.permute(1,0))
    conFa=torch.conj(Fa.permute(1,0))
    Cs=Cs.to(device)
    power_d=power_d.to(device)
    net.train()
    net.set_to_train()
    snr_w=snr_w.to(device)
    weights=weights.to(device)
    # z_batch= torch.empty(0,40,64,device=device)
    # do a training epoch over each mini-batch x returned
    # by the data loader
    for i, data in enumerate(trainloader, 0):
        # if on GPU put mini-batch into CUDA memory
        H_freq = data[0].to(device, dtype=torch.cfloat)
        H_delay= data[1].to(device, dtype=torch.cfloat)
        # do ELBO gradient and accumulate loss
        optimizer.zero_grad()

        # H_vbi=torch.zeros_like(H_delay)
        batch_len=H_freq.shape[0]

        var_idx=torch.randint(snr_w.shape[0], (batch_len,))
        var=-snr_w[var_idx]/torch.tensor(10,device=device)
        noise_var=torch.float_power(torch.tensor(10,device=device), var)   
        noise_var=noise_var.reshape(H_freq.shape[0],1,1).expand(H_freq.shape)
        noise=torch.randn(H_freq.shape,dtype=torch.complex64,device=device) 
        noise=noise*torch.sqrt(noise_var)
        H_LS=H_freq+noise
        H_LS=H_LS.to(device, dtype=torch.complex64)
        H_delay_LS=torch.matmul(conFm,H_LS) # 1000 15 64
        # H_vbi=LMMSE_fast(H_delay_LS,Filter_w,var_idx)
        H_vbi=H_delay_LS

        weight=weights[var_idx]
        
        H_vbi=torch.matmul(H_vbi,conFa)
        err=torch.square( torch.abs(H_vbi - H_delay)).mean(dim=(1,2))

        outputs, z, b = net(H_vbi,err)
        # loss = criterion(outputs, obs, prior_cost)
        loss=criterion(outputs,H_delay,z, b,weight,device=device)

        # z_batch=torch.cat((z,z_batch),dim=0)

        loss.backward(retain_graph=True)
        optimizer.step()

        # print statistics
        running_loss += loss.item() 
    # with torch.no_grad():
    #     Cs=update_Cs(z_batch)

    # return epoch loss
    normalizer_train = len(trainloader)
    total_epoch_loss_train = running_loss / normalizer_train
    return total_epoch_loss_train



def test_loop(model,dataloader,  loss_fn, Fm,Fa, device = torch.device('cpu')):
    Fm=Fm.to(device)
    Fa=Fa.to(device)
    num_batches = len(dataloader)
    test_loss = 0
    model.eval()
    pr_w = 0
    model.set_to_inference()

    with torch.no_grad():
        for i, data in enumerate(dataloader, 0):
            # get the inputs; data is a list of [inputs, labels]
            inputs, labels = data[0].to(device, dtype=torch.cfloat), data[1].to(device, dtype=torch.cfloat)
            # outputs = model.reconstruct_img(inputs)
            err=torch.square( torch.abs(inputs - labels)).mean(dim=(1,2))
            outputs, z, b =model(inputs,err)
            
            # outputs=post_NN_module(outputs,Fm,Fa)
            loss = loss_fn(outputs, labels)
            test_loss += loss.item()
            pr_w = pr_w+ torch.mean(b)

    # test_loss = test_loss/size/M/K
    test_loss = test_loss/num_batches
    pr_w = pr_w/num_batches
    
    return test_loss, pr_w

def train_criterion(output, obs, z, b,weight,device = torch.device('cpu')):
    # loss = torch.square(output - target).sum()
    output=torch.view_as_real(output)  # 1000 96 64 2
    obs_sp=torch.view_as_real(obs)
    mse1 = torch.square(output - obs_sp).sum(dim=(1,2,3))
    mse2 = torch.square( obs_sp).sum(dim=(1,2,3))
    mse = (torch.div(mse1,mse2)*weight).mean()
    
    reg1 = torch.abs(z) # 1000 15 32 
    b1 = b+torch.tensor(1e-10,device=device)
    reg1 = torch.div(reg1,b1).mean()

    Pz = torch.square(torch.abs(z)).mean()
    Pz = torch.tensor(1.0)/Pz
    reg2 = Pz*10-9*torch.log(Pz)

    # pr1=b.mean()
    # pr0=torch.tensor(1,device=device)-pr1
    # pri_1=torch.tensor(1/10,device=device)
    # pri_0=torch.tensor(1-1/10,device=device)
    # pr_logi=pr1/pr0*pri_0/pri_1
    # reg3 = torch.log(pr_logi)
    # reg3 = reg3*(reg3>0)+(pr_logi-1)*(reg3<=0)/torch.tensor(10,device=device)

    p_x=torch.square(torch.abs(obs))
    Pobs = torch.square(obs_sp).sum(dim=(3))
    x_thr=Pobs.mean(dim=(1,2))/torch.tensor(100,device=device)
    x_thr2=x_thr.reshape(x_thr.shape[0],1,1).expand(-1,p_x.shape[1],p_x.shape[2])
    sp_x=((p_x>x_thr2)*1.0).mean(dim=[1,2])

    pr1=b.mean(dim=(1,2))
    pr0=torch.tensor(1,device=device)-pr1
    pri_1=sp_x
    pri_0=torch.tensor(1,device=device)-pri_1
    pr_logi=pr1/pr0*pri_0/pri_1
    reg3 = torch.log(pr_logi)
    reg3 = reg3*(reg3>0)+(pr_logi-1)*(reg3<=0)/torch.tensor(10,device=device)
    reg3 = reg3.mean()

    loss_reg = reg1+reg2+reg3

    loss = 200*mse+loss_reg
    return loss

def test_criterion(output, target):
    # loss = torch.square(output - target).sum()
    output=torch.view_as_real(output)
    target=torch.view_as_real(target)
    loss1 = torch.square(output - target).sum(dim=(1,2,3))
    loss2 = torch.square( target).sum(dim=(1,2,3))
    loss = torch.div(loss1,loss2).mean()           
    return loss

if __name__ == '__main__':
    device = torch.device('cuda:1' if torch.cuda.is_available() else 'cpu')
    # device = 'cpu'
    print(device)
    print('Using device:', device)

    # PATH_train='./no_GB/save_result/dmrs_channel_pre_cnn_1e5.pt'
    # PATH_train='./major/save_data/dmrs_channel_3gpp_pre_cnn_1e5.pt'
    PATH_train= './major/save_data/major_channel_cdl_off_rddelaypre_cnn_1e5.pt'
    H_home=torch.load(PATH_train)['H_home'] # 1000 128 32

    Fm=torch.load(PATH_train)['Fm']
    conFm=torch.conj(Fm.permute(1,0))
    Fa=torch.load(PATH_train)['Fa']
    # Fa=torch.kron(torch.eye(2),Fa.contiguous())
    conFa=torch.conj(Fa.permute(1,0))
    # Cs=torch.load(PATH_train)['Cs']
    # power_d=torch.load(PATH_train)['power_d']

    H_freq=H_home # 1000 6 64
    H_size=H_freq.shape  # 1000 96 64 

    H_delay=torch.matmul(conFm,H_home) # 1000 36 32
    H_delay=torch.matmul(H_delay,conFa) # 1000 36 32

    batch_size_train = 256
    train_dataset =CustomImageDataset(H_input=H_freq, H_label=H_delay)
    trainloader = DataLoader(dataset=train_dataset, batch_size=batch_size_train, 
                                    shuffle=True, num_workers=4, pin_memory=True,drop_last=True)


    
    # H_test=torch.load(PATH_train)['H_test']
    PATH_t = './major/matlabcode/dmrs_channel_test_CDLA_sub60_1e3.mat'
    H_test=mat73.loadmat(PATH_t)['H_home_batch']
    H_test=torch.from_numpy(H_test)
    H_test=H_test.permute(2,0,1); # 1000 96 64
    H_test_w=H_test.cfloat()
    H_test=H_test_w[0:2000,:,:]


    Ht_delay=torch.matmul(conFm,H_test) # 1000 36 32
    Ht_delay=torch.matmul(Ht_delay,conFa) # 1000 36 32

    PATH_tpara = './major/matlabcode/major_online_vbitrain_CDLA_sub60_test_1e5.mat'
    Cs = mat73.loadmat(PATH_tpara)['cov_H_mix']  # load from .mat
    Cs=torch.from_numpy(Cs).cfloat()
    power_d = mat73.loadmat(PATH_tpara)['power_d']  # load from .mat
    power_d =torch.from_numpy(power_d).cfloat()


    PATH_filter = './major/save_result/LMSE_filter_rddlay_snr_n20_5.pt'
    snr_w=torch.load(PATH_filter)['snr_w'] 
    Filter_w=torch.load(PATH_filter)['Filter_w']  
    weights=torch.load(PATH_filter)['weights'] 
    # R_inv= torch.load(PATH_filter)['R_inv'] 

    snr_t_idx=torch.round(torch.linspace(0,snr_w.shape[0]-1, 10)).to(torch.int64)
    snr_t=snr_w[snr_t_idx]
    print(snr_t)
    weight_test=weights[snr_t_idx]
    

    test_snr=snr_t[0]
    Ht_vbi=pre_module(H_test,test_snr,Fm,Cs,power_d)
    Ht_est = torch.matmul(Fm,Ht_vbi)
    err0=NMSE(Ht_est, H_test)
    print(err0)
    Ht_vbi=torch.matmul(Ht_vbi,conFa)  
    batch_size_test = 1000
    test_dataset =CustomImageDataset(H_input=Ht_vbi, H_label=Ht_delay)
    testloader_0 = DataLoader(dataset=test_dataset, batch_size=batch_size_test, 
                                    shuffle=False, num_workers=0, pin_memory=True,drop_last=True)

    test_snr=snr_t[1]
    Ht_vbi=pre_module(H_test,test_snr,Fm,Cs,power_d)
    Ht_est = torch.matmul(Fm,Ht_vbi)
    err1=NMSE(Ht_est, H_test)
    print(err1)
    Ht_vbi=torch.matmul(Ht_vbi,conFa)  
    batch_size_test = 1000
    test_dataset =CustomImageDataset(H_input=Ht_vbi, H_label=Ht_delay)
    testloader_1 = DataLoader(dataset=test_dataset, batch_size=batch_size_test, 
                                    shuffle=False, num_workers=0, pin_memory=True,drop_last=True)                                

    test_snr=snr_t[2]
    Ht_vbi=pre_module(H_test,test_snr,Fm,Cs,power_d)
    Ht_est = torch.matmul(Fm,Ht_vbi)
    err2=NMSE(Ht_est, H_test)
    print(err2)
    Ht_vbi=torch.matmul(Ht_vbi,conFa)  
    batch_size_test = 1000
    test_dataset =CustomImageDataset(H_input=Ht_vbi, H_label=Ht_delay)
    testloader_2 = DataLoader(dataset=test_dataset, batch_size=batch_size_test, 
                                    shuffle=False, num_workers=0, pin_memory=True,drop_last=True)

    test_snr=snr_t[3]
    Ht_vbi=pre_module(H_test,test_snr,Fm,Cs,power_d)
    Ht_est = torch.matmul(Fm,Ht_vbi)
    err3=NMSE(Ht_est, H_test)
    print(err3)
    Ht_vbi=torch.matmul(Ht_vbi,conFa)  
    batch_size_test = 1000
    test_dataset =CustomImageDataset(H_input=Ht_vbi, H_label=Ht_delay)
    testloader_3 = DataLoader(dataset=test_dataset, batch_size=batch_size_test, 
                                    shuffle=False, num_workers=0, pin_memory=True,drop_last=True)                                
    
    test_snr=snr_t[4]
    Ht_vbi=pre_module(H_test,test_snr,Fm,Cs,power_d)
    Ht_est = torch.matmul(Fm,Ht_vbi)
    err4=NMSE(Ht_est, H_test)
    print(err4)
    Ht_vbi=torch.matmul(Ht_vbi,conFa)  
    batch_size_test = 1000
    test_dataset =CustomImageDataset(H_input=Ht_vbi, H_label=Ht_delay)
    testloader_4 = DataLoader(dataset=test_dataset, batch_size=batch_size_test, 
                                    shuffle=False, num_workers=0, pin_memory=True,drop_last=True)


    test_snr=snr_t[5]
    Ht_vbi=pre_module(H_test,test_snr,Fm,Cs,power_d)
    Ht_est = torch.matmul(Fm,Ht_vbi)
    err5=NMSE(Ht_est, H_test)
    print(err5)
    Ht_vbi=torch.matmul(Ht_vbi,conFa)  
    batch_size_test = 1000
    test_dataset =CustomImageDataset(H_input=Ht_vbi, H_label=Ht_delay)
    testloader_5 = DataLoader(dataset=test_dataset, batch_size=batch_size_test, 
                                    shuffle=False, num_workers=0, pin_memory=True,drop_last=True)                                
    
    test_snr=snr_t[6]
    Ht_vbi=pre_module(H_test,test_snr,Fm,Cs,power_d)
    Ht_est = torch.matmul(Fm,Ht_vbi)
    err6=NMSE(Ht_est, H_test)
    print(err6)
    Ht_vbi=torch.matmul(Ht_vbi,conFa)  
    batch_size_test = 1000
    test_dataset =CustomImageDataset(H_input=Ht_vbi, H_label=Ht_delay)
    testloader_6 = DataLoader(dataset=test_dataset, batch_size=batch_size_test, 
                                    shuffle=False, num_workers=0, pin_memory=True,drop_last=True)
    
    test_snr=snr_t[7]
    Ht_vbi=pre_module(H_test,test_snr,Fm,Cs,power_d)
    Ht_est = torch.matmul(Fm,Ht_vbi)
    err7=NMSE(Ht_est, H_test)
    print(err7)
    Ht_vbi=torch.matmul(Ht_vbi,conFa)  
    batch_size_test = 1000
    test_dataset =CustomImageDataset(H_input=Ht_vbi, H_label=Ht_delay)
    testloader_7 = DataLoader(dataset=test_dataset, batch_size=batch_size_test, 
                                    shuffle=False, num_workers=0, pin_memory=True,drop_last=True)
    test_snr=snr_t[8]
    Ht_vbi=pre_module(H_test,test_snr,Fm,Cs,power_d)
    Ht_est = torch.matmul(Fm,Ht_vbi)
    err8=NMSE(Ht_est, H_test)
    print(err8)
    Ht_vbi=torch.matmul(Ht_vbi,conFa)  
    batch_size_test = 1000
    test_dataset =CustomImageDataset(H_input=Ht_vbi, H_label=Ht_delay)
    testloader_8 = DataLoader(dataset=test_dataset, batch_size=batch_size_test, 
                                    shuffle=False, num_workers=0, pin_memory=True,drop_last=True)
    test_snr=snr_t[9]
    Ht_vbi=pre_module(H_test,test_snr,Fm,Cs,power_d)
    Ht_est = torch.matmul(Fm,Ht_vbi)
    err9=NMSE(Ht_est, H_test)
    print(err9)
    Ht_vbi=torch.matmul(Ht_vbi,conFa)  
    batch_size_test = 1000
    test_dataset =CustomImageDataset(H_input=Ht_vbi, H_label=Ht_delay)
    testloader_9 = DataLoader(dataset=test_dataset, batch_size=batch_size_test, 
                                    shuffle=False, num_workers=0, pin_memory=True,drop_last=True)                                

    err=torch.tensor([err0, err1, err2, err3,err4, err5, err6, err7,err8, err9])   
    print(err*weight_test)

    Fm=Fm.to(device)
    Fa=Fa.to(device)
    vae = VAE(device=device)
    vae.to(device)

    

    learning_rate = 1e-4
    NUM_EPOCHS = 20000
    # criterion = nn.MSELoss(reduction='sum')
    optimizer = optim.Adam(vae.parameters(), lr=learning_rate, betas=(0.9, 0.999), 
                                eps=1e-08, weight_decay=0, amsgrad=False)

 
    # Model_PATH = "./major/robust_noici_CNN_model_rddelay_novbi_para_n20_10_cdlA_testn20_0.pt"
    Model_PATH = "./major/robust_noici_CNN_model_rddelay_novbi_para_n20_3_cdlA.pt"

    PATH_tpara = './major/matlabcode/offline_vbitrain_major_sub60_1e5.mat'
    # PATH_tpara = './major/matlabcode/major_online_vbitrain_CDLB_sub60_test_1e3.mat'
    FCs = mat73.loadmat(PATH_tpara)['cov_H_mix']  # load from .mat
    FCs=torch.from_numpy(FCs).cfloat()
    Fpower_d = mat73.loadmat(PATH_tpara)['power_d']  # load from .mat
    Fpower_d =torch.from_numpy(Fpower_d).cfloat()


    TEST_FREQUENCY = 1
    min_valid_loss = np.inf
    # the_last_loss = 0.0
    patience = 1000
    trigger_times = 0
    # training loop
    # Cs = torch.eye(64,device=device,dtype=torch.cfloat)
    for epoch in range(NUM_EPOCHS):
        total_epoch_loss_train = train(vae,trainloader, train_criterion,optimizer, Fm,Fa,FCs,Fpower_d,
                                snr_w,weights,device=device)
        
        print("[epoch %03d]  average training loss: %.4f" % (epoch, total_epoch_loss_train))

        if epoch % TEST_FREQUENCY == 0:
            # report test diagnostics
            tloss_0, pr_w0=test_loop(vae, testloader_0, test_criterion,Fm,Fa, device=device)
            tloss_1, pr_w=test_loop(vae, testloader_1, test_criterion,Fm,Fa, device=device)
            tloss_2, pr_w=test_loop(vae, testloader_2, test_criterion,Fm,Fa, device=device)
            tloss_3, pr_w=test_loop(vae, testloader_3, test_criterion,Fm,Fa, device=device)
            tloss_4, pr_w=test_loop(vae, testloader_4, test_criterion,Fm,Fa, device=device)
            tloss_5, pr_w=test_loop(vae, testloader_5, test_criterion,Fm,Fa, device=device)
            tloss_6, pr_w=test_loop(vae, testloader_6, test_criterion,Fm,Fa, device=device)
            tloss_7, pr_w=test_loop(vae, testloader_7, test_criterion,Fm,Fa, device=device)
            tloss_8, pr_w=test_loop(vae, testloader_8, test_criterion,Fm,Fa, device=device)
            tloss_9, pr_w9=test_loop(vae, testloader_9, test_criterion,Fm,Fa, device=device)

            # total_epoch_loss_test = evaluate(svi, test_loader, use_cuda=USE_CUDA)

            print("[epoch %03d] snr=%.2f, test loss: %.4f" % (epoch, snr_t[0], tloss_0/err0))
            print("[epoch %03d] snr=%.2f, test loss: %.4f" % (epoch, snr_t[1], tloss_1/err1))
            print("[epoch %03d] snr=%.2f, test loss: %.4f" % (epoch, snr_t[2], tloss_2/err2))
            print("[epoch %03d] snr=%.2f, test loss: %.4f" % (epoch, snr_t[3], tloss_3/err3))
            print("[epoch %03d] snr=%.2f, test loss: %.4f" % (epoch, snr_t[4], tloss_4/err4))
            print("[epoch %03d] snr=%.2f, test loss: %.4f" % (epoch, snr_t[5], tloss_5/err5))
            print("[epoch %03d] snr=%.2f, test loss: %.4f" % (epoch, snr_t[6], tloss_6/err6))
            print("[epoch %03d] snr=%.2f, test loss: %.4f" % (epoch, snr_t[7], tloss_7/err7))
            print("[epoch %03d] snr=%.2f, test loss: %.4f" % (epoch, snr_t[8], tloss_8/err8))
            print("[epoch %03d] snr=%.2f, test loss: %.4f" % (epoch, snr_t[9], tloss_9/err9))
            print("[epoch %03d] sparsity level: %.4f" % (epoch, pr_w))

            # test_loss=torch.max(torch.tensor([tloss_0, tloss_1, tloss_2,tloss_3,
            #     tloss_4, tloss_5, tloss_6,tloss_7,tloss_8, tloss_9])*weight_test)
            # test_loss=torch.mean(torch.tensor([tloss_0/err0, tloss_1/err1, tloss_2/err2,tloss_3/err3,
            #     tloss_4/err4, tloss_5/err5, tloss_6/err6,tloss_7/err7,tloss_8/err8, tloss_9/err9]))
            test_loss=torch.mean(torch.tensor([tloss_0/err0, tloss_1/err1, tloss_2/err2,tloss_3/err3,
                tloss_4/err4, tloss_5/err5, tloss_6/err6]))
            # test_loss=torch.mean(torch.tensor([tloss_3/err3,
            #     tloss_4/err4, tloss_5/err5, tloss_6/err6,tloss_7/err7]))
            # test_loss=torch.mean(torch.tensor([tloss_1/err1, tloss_2/err2,tloss_3/err3,
            #     tloss_4/err4, tloss_5/err5, tloss_6/err6]))
            print("[epoch %03d] test ratio: %.4f" % (epoch, test_loss))
            print("[epoch %03d] sparsity level 0: %.4f" % (epoch, pr_w0))
            print("[epoch %03d] sparsity level 9: %.4f" % (epoch, pr_w9))

        if test_loss-min_valid_loss > 0.0001:
            trigger_times += 1
            # print('trigger times:', trigger_times)
            print("trigger_times: %.1f" % trigger_times)
            if trigger_times >= patience:
                print('Early stopping due to validation loss increase!\n')
                break
        else:
            # print('trigger times: 0')
            trigger_times = 0
        

        if min_valid_loss > test_loss:
            min_valid_loss = test_loss
            # Save
            torch.save(vae.state_dict(), Model_PATH)
            print("model updated")
   