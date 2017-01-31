function [ ux ] = essfm( Pavg,ux,C,ch )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

%         
        ux = sig.FIELDX;
        dz = ch.Lf;
        
        if abs(ch.alphalin*dz) > 1e-6
            Leff     = -(1-exp(ch.alphalin*dz))/ch.alphalin;
        else
            Leff     = dz;
        end
        
        z    = dz*(0:1-1);
        
        xi   = -ch.gamma*Leff*exp(-ch.alphalin*z)*Pavg;
        
        halfdz      = dz/2;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %                         HALF DZ GVD                        %
        %                              DZ SPM                        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        ux = lin_step(-ch.beta*halfdz,ux);
        ux = essfm_nl_step(-xi(1)*C,ux);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        for i=2:1
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %                           DZ GVD                       %
            %                           DZ SPM                       %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            ux = lin_step(-ch.beta*dz,ux);
            ux = essfm_nl_step(-xi(i)*C,ux);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %                      LAST HALF DZ GVD                      %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        ux = lin_step(-ch.beta*halfdz,ux);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
         sig.FIELDX = ux;

end

 function ux         = lin_step(betaxdz,ux)
        
        Hf = exp(-1i*betaxdz);
        ux = ifft( fft(ux) .* Hf);
        
    end

    function [ux]       = essfm_nl_step(C,ux)
        
        pow = real(ux).^2 + imag(ux).^2 ;
        
        M=length(pow);
        N=length(C);
        per_pow =[pow(M-N+2:M);pow;pow(1:N-1)];
        
        CC=[flipud(C);C(2:end)];
        theta=conv(per_pow,CC,'valid');
        
        ux = ux.*exp(1i*theta);
        
    end
