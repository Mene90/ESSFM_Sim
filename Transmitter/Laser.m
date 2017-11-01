classdef Laser
    %LASERSOURCE Summary of this class goes here
    %   Detailed explanation goes here
    
       
    methods (Static = true)
        function E = GetLaserSource(Pin,sig,lam,spac)
            CLIGHT = 299792458;
            E   = ones(sig.NSYMB*sig.NT,1); %*sqrt(Pin);
            La = length(lam);               % init wavelengths
            if (La == 1) && (sig.NCH > 1)
                if ~exist('spac','var')
                    error('missing the channel-spacing SPAC');
                end
                lamt = zeros(1,sig.NCH);
                if (spac == 0)
                    if(mod(sig.NCH,2))
                        deltafn = (0:sig.SYMBOLRATE:sig.SYMBOLRATE*(sig.NCH-1))-((sig.SYMBOLRATE/2)*(sig.NCH-1));
                    else
                        deltafn = (0:sig.SYMBOLRATE:sig.SYMBOLRATE*(sig.NCH-1))-((sig.SYMBOLRATE/2)*(sig.NCH-1));
                    end
                    lamt    =  -(deltafn./CLIGHT-1/lam).^-1;
                else
                    for chan = 1:sig.NCH
                        lamt(chan)=lam+spac*(chan-(sig.NCH+1)/2);
                    end
                end
            elseif La ~= sig.NCH
                error('wrong length for LAM (must be 1 or # channels)');
            else
                lamt = lam;
            end
            
            
            sig.LAMBDA = lamt;  % global and unique definition of wavelengths
            
        end
        
        % Add gaussian complex white noise (ASE)?
    end
    
end

