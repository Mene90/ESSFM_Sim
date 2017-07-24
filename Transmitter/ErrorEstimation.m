classdef ErrorEstimation
    %BER Summary of this class goes here
    %   Detailed explanation goes here
        
     methods (Static = true)
         
         function ber = BER(sig)
             
             patmatx_tx = samp2pat(angle(sig.FIELDX_TX));
             patmatx    = samp2pat(angle(sig.FIELDX));
             if any(sig.FIELDY)
                patmaty_tx = samp2pat(angle(sig.FIELDY_TX(1:sig.NT:end)));     
                patmaty    = samp2pat(angle(sig.FIELDY(1:sig.NT:end)));
                ber        = (ber_aux(patmatx,patmatx_tx)  +...
                              ber_aux(patmaty,patmaty_tx)) * 0.5;
             else
                ber        =  ber_aux(patmatx,patmatx_tx);
             end
             
         end
         
         function mse = MSE(sig)
             if any(sig.FIELDY)
                 rx = [sig.FIELDX'; sig.FIELDY'];
                 tx = [sig.FIELDX_TX'; sig.FIELDY_TX'];
                 mse = mse_aux(rx,tx);
             else
                 rx = sig.FIELDX;
                 tx = sig.FIELDX_TX;
                 mse = mean(abs(rx-tx).^2);
             end 
         end
    end
    
end

function [ avg_ber ] = ber_aux(patmat_hat,patmat_rx )
    avg_ber = (mean(patmat_hat(:,1)~=patmat_rx(:,1))  +...
               mean(patmat_hat(:,2)~=patmat_rx(:,2))) * 0.5;
end

function [ avg_mse ] = mse_aux(rx,tx)
       sum  = 0;
       for i = size(rx,2)   
          sum = sum + abs(rx(:,i)-tx(:,i)).^2;
       end
       avg_mse = sum/i;
end

