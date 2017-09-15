function [ irate ] = AIR(sig, modtype)
        xk = sig.FIELDX_TX;
        yk = sig.FIELDX;
        
        sigmaxy  = sum(yk.*conj(xk));
        sigmax2  = sum(abs(xk).^2);
        sigmay2  = sum(abs(yk).^2);
        
        h       = sigmaxy./sigmax2;
        
        zk      = yk./h;
        N       = length(xk);
        
        SNR     = abs(sigmaxy)^2/(sigmax2*sigmay2 - abs(sigmaxy)^2);
        
        if(strcmp(modtype,'Gaussian'))
        
            irate = log2(SNR + 1) - SNR/(N*log(2)) * sum(abs(zk-xk).^2 - abs(zk).^2/(SNR + 1)); 

        else
             error('Only gaussian modulation is implemented');
        end
end

