classdef Ampliflat
    
    properties 
        HPLANCK = 6.62606896e-34;       % Planck's constant [J*s]
        CLIGHT = 299792458;             % speed of light [m/s]
    end
    
    properties 
        N0;
    end
    
    methods
        
        function amp = Ampliflat(Pavg,ch,Gm1,etasp,type,Nspan)
          
          if (strcmp(type,'Raman'))            
            ktusual     = etasp;
            amp.N0      = ch.Lf*ch.alphalin*amp.HPLANCK*amp.CLIGHT/(ch.lambda * 1e-9)*ktusual/Pavg*Nspan;
            ch.alphadB  = 0;
            ch.alphalin = 0;            
          else                                                 
            amp.N0  = etasp*(Gm1)*amp.HPLANCK*amp.CLIGHT...
                        /(ch.lambda * 1e-9)*Nspan;
          end
          
        end
                               
        function AddNoise(amp,sig)
            
            nfft  = sig.NSYMB * sig.NT;
            Df    = sig.SYMBOLRATE*sig.NT*1e+9;
            sigma = ones(sig.NSYMB * sig.NT,1)*sqrt(0.5*amp.N0*Df);
                        
            noiseX = sigma.* (randn(nfft,1)+1i*randn(nfft,1));
            set(sig,'FIELDX', sig.FIELDX + noiseX);
            
            if(any(sig.FIELDY))
                noiseY = sigma.* (randn(nfft,1)+1i*randn(nfft,1));
                set(sig,'FIELDY', sig.FIELDY + noiseY);
            end
            
        end
        
        function AddNoiseXspan(amp,sig,Nspan)
            
            nfft  = sig.NSYMB * sig.NT;
            Df    = sig.SYMBOLRATE*sig.NT*1e+9;
            sigma = ones(sig.NSYMB * sig.NT,1)*sqrt(0.5*amp.N0/Nspan*Df);
                        
            noiseX = sigma.* (randn(nfft,1)+1i*randn(nfft,1));
            set(sig,'FIELDX', sig.FIELDX + noiseX);
            
            if(any(sig.FIELDY))
                noiseY = sigma.* (randn(nfft,1)+1i*randn(nfft,1));
                set(sig,'FIELDY', sig.FIELDY + noiseY);
            end
            
        end
        
    end
    
end

