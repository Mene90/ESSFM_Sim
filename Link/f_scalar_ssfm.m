function [ sig ] = f_scalar_ssfm( ch,Pavg,sig )
    %SCALAR_SSFM Summary of this function goes here
    %   Detailed explanation goes here

    omega = 2*pi*sig.SYMBOLRATE*sig.FN*1e9;            % [rad/s]
    beta  = 0.5*omega.^2*ch.b2 + omega.^3*ch.b3/6;

    ux    = sig.FIELDX;%;get(sig,'FIELDX');

    if abs(ch.alphalin*ch.dz) > 1e-6
        Leff     = (1-exp(-ch.alphalin*ch.dz))/ch.alphalin;
    else
        Leff     = ch.dz;
    end
    steps = (0:ch.nstep-1);
    z    = ch.dz * steps;

    xi   = ch.gamma*Leff*exp(-ch.alphalin*z)*Pavg;

    halfdz      = ch.dz/2;

    
    %                           HALF DZ GVD                  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    %
    %                                DZ SPM                      %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ux = scalar_lin_step(beta*halfdz,ux);
    ux = scalar_nl_step(ux,xi(1));
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    for i=2:ch.nstep

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %                           DZ GVD                       %
        %                           DZ SPM                       %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        ux = scalar_lin_step(beta*ch.dz,ux);
        ux = scalar_nl_step(ux,xi(i));
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %                      LAST HALF DZ GVD                      %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ux = scalar_lin_step(beta*halfdz,ux);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    sig.FIELDX = ux;%set(sig,'FIELDX',ux);


end

function ux  = scalar_lin_step(betaxdz,ux)

    Hf = exp(-1i*betaxdz);
    ux = ifft( fft(ux) .* Hf);

end

function ux  = scalar_nl_step(ux,xi)

    pow = real(ux).^2 + imag(ux).^2;
    ux  = ux.*exp(-1i*xi.*pow);

end