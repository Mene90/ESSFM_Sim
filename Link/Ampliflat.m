classdef Ampliflat
    %AMPLIFLAT Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        HPLANCK = 6.62606896e-34;       % Planck's constant [J*s]
        CLIGHT = 299792458;             % speed of light [m/s]
    end
    
    properties (SetAccess = private)
        N0;
    end
    
    methods
        
        function amp = Ampliflat(Pavg,ch,G,etasp)
            Gm1 =(10^(G*0.1)-1.d0);                                    % (G-1)
            amp.N0  = etasp*(Gm1)*amp.HPLANCK*amp.CLIGHT/(ch.lambda * 1e-9)/Pavg;
            
%             N0=Gm1*obj.HPLANCK*obj.CLIGHT/(ch.lambda * 1e-9)*etasp;  % Accumulate ASE Noise PSD
% %             SNR=Pavg/(ch.signal.SYMBOLRATE * 1e+9)/N0;
%             SNR=1/(ch.signal.SYMBOLRATE * 1e+9)/N0;
%             obj.sigma=sqrt(ch.signal.NT*0.5/SNR);
        end
        
        function sig = amplifie(amp,sig,G)
            
             if not(isa(sig,'Signal'))
                error(['sig parameter has to be an instance',...
                    ' of Signal class']);
             end
            
             sig.FIELDX = sig.FIELDX*sqrt(10^(G*0.1));
             sig.FIELDY = sig.FIELDY*sqrt(10^(G*0.1));
        end
                
        function sig = AddNoise(amp,sig)
            
            nfft = sig.NSYMB * sig.NT;
            
            Df    = sig.SYMBOLRATE*sig.NT*1e+9;
            sigma = ones(nfft,1)*sqrt(0.5*amp.N0*Df);
                        
            noiseX = sigma.* (randn(nfft,1)+1i*randn(nfft,1));
            sig.FIELDX = sig.FIELDX + noiseX;
            
            if(not(isempty(sig.FIELDY)))
                noiseY = sigma.* (randn(nfft,1)+1i*randn(nfft,1));
                sig.FIELDY = sig.FIELDY + noiseY;
            end
        end
        
        function [ux,uy] = gpu_AddNoise(amp,sig,ux,uy)
            
            nfft = sig.NSYMB * sig.NT;
            
            Df    = sig.SYMBOLRATE*sig.NT*1e+9;
            sigma = ones(nfft,1)*sqrt(0.5*amp.N0*Df);
                        
            noiseX = sigma.* (randn(nfft,1)+1i*randn(nfft,1));
            ux = ux + noiseX;
            
            if(not(isempty(sig.FIELDY)))
                noiseY = sigma.* (randn(nfft,1)+1i*randn(nfft,1));
                uy = uy + noiseY;
            end
        end
        
    end
    
end

