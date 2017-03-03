function [ max_snr ] = SSFM_MAX_SNR( Nstep,sym_length,n_prop_steps,R,etasp,Nspan)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                         Link parameters                                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
LL        = 1e5;                % length [m]
alphadB   = 0.2;                  % attenuation [dB/km]
aeff      = 85;                   % effective area [um^2]
n2        = 2.5e-20;              % nonlinear index [m^2/W]
lambda    = 1550;                 % wavelength [nm] @ dispersion
D         = 17;                   % dispersion [ps/nm/km] @ wavelength
S         = 0;                    % slope [ps/nm^2/km] @ wavelength

Ns_prop   = n_prop_steps;         % number of SSFM propagation step
pmd       = false;                % pmd enable/disable

ch        = Channel(LL,alphadB,lambda,aeff,n2,D,S,Ns_prop,pmd);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                         DSP parameters                                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Ns_bprop  = Nstep;                  % SSFM and ESSFM backpropagation steps
dsp       = DSP(ch,Ns_bprop);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                      Global Signal parameters                          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
symbrate  = R;                  % symbol rate [Gbaud]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                         Signal parameters                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Nsymb     = sym_length;                % number of symbols
Nt        = 2;                         % points x symbol
sig       = Signal(Nsymb,Nt,symbrate);

max_snr   = -ssfm_maxSNR(ch,dsp,sig,etasp,Nspan);
end

