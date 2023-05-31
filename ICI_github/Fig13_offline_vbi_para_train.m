close all
clear all

% load('dmrs_channel_home_mixmodel_1e4.mat','H_home_batch');
% load('dmrs_channel_home_3gppmodel_1e5.mat','H_home_batch');
% load('dmrs_channel_256_home_mixmodel_1e4.mat','H_home_batch');
% load('dmrs_channel_home_mixmodel_major_1e5.mat','H_home_batch');
load('dmrs_channel_home_mixmodel_major_sub60_1e5.mat','H_home_batch');
[Nd,M,ite]=size(H_home_batch);
% Fa=kron((fftshift(dftmtx(8),1)./sqrt(8)).',fftshift(dftmtx(4),2)./sqrt(4)); %% Angular basis
Fa=kron((fftshift(dftmtx(8),1)./sqrt(8)),fftshift(dftmtx(4),1)./sqrt(4));
Lc=round(144/1024*Nd);
Fm=dftmtx(Nd)/sqrt(Nd);
Fm(:,(Lc+1):(end-Lc))=[];
Fm = circshift(Fm,Lc,2); %% Basis for delay-domain signal
L=2*Lc;
Fi=dftmtx(Nd)/sqrt(Nd);
Fi=Fi(:,(Lc+1):(end-Lc)); %% Basis for delay-domain interference
%%
H=reshape(H_home_batch,[Nd,M*ite]);
H_delay=Fm'*H;
% P_delay=mean(abs(H_delay).^2,2);
% figure
% semilogy(1:L,P_delay,'r.');
% grid on
H_delay=reshape(H_delay,[L,M,ite]);
H_delay=permute(H_delay,[2,3,1]); %% M,ite,L
%%
R_w=zeros(M,M,L);
for l0=1:L
    H_ell=H_delay(:,:,l0);
    R_tmp=zeros(M,M);
    for i0=1:ite
        H_tmp= H_ell(:,i0);
        R_tmp=R_tmp+H_tmp*H_tmp';
    end
    R_w(:,:,l0)=R_tmp;
end
%%
Eta=eye(M);
alpha=M*ite;
loop=1;
while(1)
    beta_w=zeros(1,L);
    for l0=1:L
        beta_w(l0)=trace(R_w(:,:,l0)*Eta);
    end
    Eta_inv=zeros(M,M);
    for l0=1:L
        Eta_inv=Eta_inv+alpha/beta_w(l0)*R_w(:,:,l0);
    end
    Eta_inv=Eta_inv./L./ite;
    %%
    ebs=trace(Eta_inv)/M;
    Eta_inv=Eta_inv./ebs;
    beta_w=beta_w.*ebs;
    %%
    Eta_new=inv(Eta_inv);
    flag=norm(Eta_new-Eta)
    if flag<1e-5
        break
    end
    Eta=Eta_new;
    loop=loop+1;
end
cov_H_mix=Eta_inv;
power_d=beta_w./alpha;
R_w=R_w./ite;
%%
%%
% save('offline_vbitrain_1e4.mat','cov_H_mix','power_d','Fm','Fa','Fi','R_w');
% save('offline_vbitrain_1e4.mat','cov_H_mix','power_d','Fm','Fa','Fi','R_w');
% save('offline_256_vbitrain_1e4.mat','cov_H_mix','power_d','Fm','Fa','Fi','R_w');
% save('offline_vbitrain_3gpp_1e5.mat','cov_H_mix','power_d','Fm','Fa','Fi','R_w');
% save('offline_vbitrain_major_1e5.mat','cov_H_mix','power_d','Fm','Fa','Fi','R_w');
save('offline_vbitrain_major_sub60_1e5.mat','cov_H_mix','power_d','Fm','Fa','Fi','R_w');
%%




