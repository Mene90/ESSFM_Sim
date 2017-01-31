function [ sig ] = ampli( Pavg,ch,G,etasp )
%AMPLI Summary of this function goes here
%   Detailed explanation goes here

        HPLANCK = 6.62606896e-34;      % Planck's constant [J*s]
        CLIGHT = 299792458;      % speed of light [m/s]
        
        sig = ch.signal;
        nfr = length(sig.FIELDX);

        %Gerbio   = ch.alphadB*ch.length/1e3
        %gain     = 10^(Gerbio*0.1);

        %amplifico
        sig.FIELDX = sig.FIELDX*sqrt(G);
        sig.FIELDY = sig.FIELDY*sqrt(G);
        
        
        %add noise
        Gm1=(G-1.d0);                                 % (G-1)
        N0=Gm1*HPLANCK*CLIGHT/(ch.lambda * 1e-9)*etasp;        %Accumulate ASE Noise PSD
        SNR=Pavg/(sig.SYMBOLRATE * 1e+9)/N0;


        sigma=sqrt(sig.NT*0.5/SNR);

        sigma_matrix = ones(nfr,1) * sigma;

        noiseX = sigma_matrix ;%.* (randn(nfr,1)+1i*randn(nfr,1));
        sig.FIELDX = sig.FIELDX + noiseX;
        
        if(sig.FIELDY ~=0)
            noiseY = sigma_matrix .* (randn(nfr,1)+1i*randn(nfr,1));
            sig.FIELDY = sig.FIELDY + noiseY;
        end


end

