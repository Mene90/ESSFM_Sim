function [ avg_ber ] = ber( patmat_hat,patmat_rx )
   

    avg_ber = 0.5*(mean(patmat_hat(:,1)~=patmat_rx(:,1))+ mean(patmat_hat(:,2)~=patmat_rx(:,2)));
      

end

