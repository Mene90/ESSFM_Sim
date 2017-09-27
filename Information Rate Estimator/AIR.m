function [ irate ] = AIR(sig, modtype)
        shift = length(sig.FIELDX_TX)*0.25;
        xk = sig.FIELDX_TX(shift:length(sig.FIELDX_TX)-shift,1);
        yk = sig.FIELDX(shift:length(sig.FIELDX)-shift,1);
        
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

