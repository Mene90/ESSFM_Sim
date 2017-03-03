function y = essfm(x,Px,Tc,alf,bt2,gm,L,Ns,nl_step,C)
%ESSFM propagates a sampled signal x through a span of fiber according to
%the nonlinear Schroedinger equation with the enhanced SSFM algorithm
%   y = essfm(x,Px,Tc,alf,bt2,gm,L,Ns,tipo,C)
%   computes the output signal y given the input signal x by using the
%   enhanced symmetric split-step Fourier method (author M. Secondini)
%   The expansion bandwidth should already account for bandwidth increase
%   during propagation
%
%   *** ARGOMENTI DELLA FUNZIONE ***
%   x:      input signal (with normalized power)
%   Px:     input signal power (the actual signal will be sqrt(Px)*x) [mW]
%   Tc:     sampling time [ps]
%   alf:    fiber attenuation [dB/km]
%   bt2:    GVD coefficient [ps^2/km]
%   gm:     nonlinear coefficient [1/(mW*km)]
%   L:      length [km]
%   Ns:     number of steps
%   nl_step:   tipo di passo nonlineare
%   C:     nonlinear coefficients [1/(mW*km)]

% Evaluation of some useful parameters
N=length(x);                    % number of samples
N2=floor(N/2);
dz=L/Ns;                        % step size
a=alf*0.230258509299405*1e-3;        % attenuation parameter (1/m)
f=1./(N*Tc)*[0:N2,-(N-N2-1):-1];% frequency vector after fft

% Dispersion vectors
phi=pi*pi*bt2*dz*f.*f;
F=exp(-1i*mod(2*phi,2*pi));       % transfer function (GVD) of step dz
Fh=exp(-1i*mod(phi,2*pi));    % transfer function (GVD) of half step dz/2

% Vector of effective nonlinearity values along z
if (abs(a*dz)>1e-6)
    geff0=gm*Px*(1.0-exp(-a*dz))/a;
else
    geff0=gm*Px*dz;
end
z=dz*(0:Ns-1);
geff=geff0*exp(-a*z);


%% First step (half disp + NL)
% linear step:
y=fft(x);
y=Fh.*y;
y=ifft(y);
% nonlinear step:
switch nl_step
    case {'filt','filtr'},
        y=transpose(nlstep_min(-geff(1)*C,transpose(y)));
    case {'quadform','quadformr'},
        y=transpose(nlstep(-geff(1)*C,transpose(y)));
end

%% Intermediate steps (whole disp + NL)
for n=2:Ns,
%   linear step:    
    y=fft(y);
    y=F.*y;
    y=ifft(y);
%   nonlinear step:
    switch nl_step
        case {'filt','filtr'},
            y=transpose(nlstep_min(-geff(n)*C,transpose(y)));
        case {'quadform','quadformr'},
            y=transpose(nlstep(-geff(n)*C,transpose(y)));
    end
end

% Final step (half disp)
y=fft(y);
y=Fh.*y;
y=ifft(y);

end

