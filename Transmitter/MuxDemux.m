classdef MuxDemux
    %MUX Summary of this class goes here
    %   Detailed explanation goes here
    
    properties(Constant)
        CLIGHT = 299792458;             % speed of light [m/s]
    end
    
    methods (Static = true)
       
        
        function [muxX,muxY] = Mux(Ex,Ey,sig)
            CLIGHT = 299792458;
            lamt = sig.LAMBDA;
            maxl = max(lamt);
            minl = min(lamt);
            fnyqmin = (CLIGHT/minl-CLIGHT/maxl)/sig.SYMBOLRATE;
            
            if (sig.NT < fnyqmin) && fnyqmin ~= 0
                reply = input('The number of samples per symbol is too small to avoid aliasing. Do you want to continue? Y/N [N]: ', 's');
                if isempty(reply)
                    reply = 'N';
                end
                if ~strcmpi(reply,'Y')
                    error('number of samples per symbol is too small');
                end
            end
            
            lamc    = round(2*maxl*minl/(maxl+minl));          % central wavelength: 1/lamc = 0.5(1/maxl+1/minl)
            deltafn = CLIGHT*(1/lamc-1./lamt);                 % absolute frequency spacing [GHz]           
            minfreq = sig.FN(2)-sig.FN(1);                     % minfreq = 1/sig.NSYMB
            ndfn = round(deltafn./sig.SYMBOLRATE/minfreq);     % spacing in points
             
%             if sig.NCH == 2
%                 ndfn = ndfn/2;
%             end
            
            zfieldx = fft(Ex);
            zfieldy = fft(Ey);
            
            for kch=1:sig.NCH    % create the unique field combining the channels
                sig.FIELDX = sig.FIELDX + fastshift(zfieldx(:,kch),-ndfn(kch));
                if any(Ey)                    
                    sig.FIELDY = sig.FIELDY + fastshift(zfieldy(:,kch),-ndfn(kch));
                end
            end
            
            sig.FIELDX = ifft(sig.FIELDX);
            
            if any(sig.FIELDY)
                sig.FIELDY = ifft(sig.FIELDY);
            end   
            
            muxX = sig.FIELDX;
            muxY = sig.FIELDY;
            
        end
        
        function [zfieldx,zfieldy] = Demux (sig,Hf,nch)
            
            CLIGHT = 299792458;
            minfreq = sig.FN(2)-sig.FN(1);
            maxl=max(sig.LAMBDA);
            minl=min(sig.LAMBDA);
            lamc = round(2*maxl*minl/(maxl+minl));                   % central wavelength
            deltafn = CLIGHT*(1/lamc-1./sig.LAMBDA);                 % frequency spacing
            ndfn = round(deltafn./sig.SYMBOLRATE/minfreq);           % spacing in points
            
%             if sig.NCH == 2
%                 ndfn = ndfn/2;
%             end
            
            sig.FIELDX = fft(sig.FIELDX);
            sig.FIELDY = fft(sig.FIELDY);
            
            if(nch == 0)
                for kch=1:sig.NCH    % create the unique field combining the channels
                    zfieldx(:,kch) = fastshift(sig.FIELDX,ndfn(kch));
                    zfieldx(:,kch) = ifft(zfieldx(:,kch) .* Hf);
                    if any(sig.FIELDY)
                        zfieldy(:,kch) = fastshift(sig.FIELDY,ndfn(kch));
                        zfieldy(:,kch) = ifft(zfieldy(:,kch) .* Hf);
                    end
                end
                sig.FIELDX = ifft(sig.FIELDX);
                sig.FIELDY = ifft(sig.FIELDY);
            else
                zfieldx(:,nch) = fastshift(sig.FIELDX,ndfn(nch));
                zfieldy(:,nch) = fastshift(sig.FIELDY,ndfn(nch));

                zfieldx(:,nch) = ifft(zfieldx(:,nch) .* Hf);            
                zfieldy(:,nch) = ifft(zfieldy(:,nch) .* Hf);
                
                sig.FIELDX = ifft(sig.FIELDX);
                sig.FIELDY = ifft(sig.FIELDY);
            end
                        
        end
    
    end
end
