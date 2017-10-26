classdef Modulator
    %MODULATOR Summary of this class goes here
    %   Detailed explanation goes here
    
    properties 
    end
    
    methods (Static = true)
        
        function E = ApplyModulation(~,samples, sig, pls)
            
            Nt     = sig.NT;
            nfft   = sig.NT * sig.NSYMB;
            H    = filt(pls,sig.FN);
%             H    = myfilter('ideal',sig.FN,1,0);
            dt   = 1/Nt;
            norm = sqrt(nfft/dt/sum(abs(H).^2));
            H    = norm*H;
            
            E=zeros(nfft,1);
            E(1:Nt:nfft)=samples;
            
            Y=fft(E);
            Y=Y.*H;
            E=ifft(Y);
        end
        
        
%         function E = ApplyModulation(~, patmat_tx, sig , pls)
%             
%             Nt     = sig.NT;
%             nfft   = sig.NT * sig.NSYMB;
%             usrad2 = 1./sqrt(2.);
%             a      = usrad2*(patmat_tx(:,1)+1i*patmat_tx(:,2));
%             H      = filt(pls,sig.FN);
%             
%             %normalizzazione per avere energia unitaria
%             dt   = 1/Nt;
%             norm = sqrt(nfft/dt/sum(abs(H).^2));
%             H    = H*norm;
% 
%             %sequenza di impulsi di tipo delta (con Nxs campioni per simbolo) modulati
%             %secondo i simboli dati
%             E=zeros(1,nfft);
%             E(1:Nt:nfft)=a;
% 
%             %shaping spettrale secondo l'impulso desiderato
%             Y=fft(E);
%             Y=Y.*H;
%             E=transpose(ifft(Y));
%         
%         end
        
%         function E = GaussianModulation(~, patmat, sig , pls)
%             Nt     = sig.NT;
%             nfft   = sig.NT * sig.NSYMB;
%             H      = filt(pls,sig.FN);
%             
%             %normalizzazione per avere energia unitaria
%             dt   = 1/Nt;
%             norm = sqrt(nfft/dt/sum(abs(H).^2));
%             H    = H*norm;
% 
%             %sequenza di impulsi di tipo delta (con Nxs campioni per simbolo) modulati
%             %secondo i simboli dati
%             E=zeros(1,nfft);
%             E(1:Nt:nfft)=patmat;
% 
%             %shaping spettrale secondo l'impulso desiderato
%             Y=fft(E);
%             Y=Y.*H;
%             E=transpose(ifft(Y));
%         end
        
    end
    
end

