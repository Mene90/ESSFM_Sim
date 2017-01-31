function Eout = mz_modulator( Ein, modsig, options )
%MZ_MODULATOR Summary of this function goes here
%   Detailed explanation goes here

        exratio     = inf;

        if exist('options','var')
            
            if isfield(options,'exratio')
                exratio = options.exratio;
            end

        end


        exr_lin = 10^(-exratio/10);
        gamma = (1-sqrt(exr_lin))/(sqrt(exr_lin)+1);

        Phi_U = (modsig *  pi/2 );
        Phi_L = (modsig * -pi/2 );

        Eout  = 1j*Ein.*(exp(1i*Phi_L)-gamma*exp(1i*Phi_U))/(1+gamma);
end

