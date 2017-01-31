classdef QI_modulator
    %IQ_MODULATOR Summary of this class goes here
    %   Detailed explanation goes here
        
    methods (Static = true)
        function E = ApplyModulation(E, modsig_i, modsig_q)
            % The parameters MODSIG_I and MODSIG_Q are the in-phase 
            % and in-quadrature electrical driving signals produced 
            % by ELECTRICSOURCE.
            
            iqratio = 0;
            exratio     = [inf inf];
            
            if exist('options','var')
                if isfield(options,'iqratio');
                    iqratio = options.iqratio;
                end
                if isfield(options,'exratio');
                    exratio = options.exratio;
                    if length(exratio) ~= 2
                        error('exratio must be a vector of length2')
                    end
                end
            end
            
            
            iqratio = 10^(iqratio/20);
            sr = iqratio/(1+iqratio);
            
            Ei = E * sr;         % In Phase Field
            Eq = E * (1-sr);     % In Quadrature Field
            
            Ei = mz_modulator(Ei,modsig_i,struct('exratio',exratio(1)));
            Eq = mz_modulator(Eq,modsig_q,struct('exratio',exratio(2)));
            
            E = sqrt(2)*(Ei+Eq*exp(1i*(pi/2)));
            
        end
    end
    
end

