classdef Laser
    %LASERSOURCE Summary of this class goes here
    %   Detailed explanation goes here
    
       
    methods (Static = true)
        function E = GetLaserSource(Pin,nfft)
            E   = ones(nfft,1);%*sqrt(Pin);
        end
        
        % Add gaussian complex white noise (ASE)?
    end
    
end

