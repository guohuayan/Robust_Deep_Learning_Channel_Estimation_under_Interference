function H = batch_channel_angular(H,Fa)
    Fa=kron(eye(2),Fa);
    [Np,num_pilot,M,Data_batch]=size(H);
    H=permute(H,[1,2,4,3]); %% Np,num_pilot,Data_batch, M
    H= reshape(H,[Np*num_pilot*Data_batch, M]);
    H=H*Fa;
    H= reshape(H,[Np,num_pilot,Data_batch, M]);
    H=permute(H,[1,2,4,3]); %% Np,num_pilot, M, Data_batch
end

