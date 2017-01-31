function [ avg_ber ] = ber( patmat_hat,patmat_rx )
%BER Summary of this function goes here
%   Detailed explanation goes here
    

    avg_ber = 0.5*(mean(patmat_hat(:,1)~=patmat_rx(:,1))+ mean(patmat_hat(:,2)~=patmat_rx(:,2)));
      

end

