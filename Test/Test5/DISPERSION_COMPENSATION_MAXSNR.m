function [ max_snr ] = DISPERSION_COMPENSATION_MAXSNR( Nstep,sym_length,n_prop_steps,R,etasp,Nspan)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                      Global Signal parameters                          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
symbrate  = R;                   % symbol rate [Gbaud]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                         Link parameters                                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
LL        = 1e5;                  % length [m]
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
%                         Amplifier parameters                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Gerbio    = ch.alphalin/(log(10)*1e-4)*ch.Lf*1e-3;
ampli.G   = Gerbio;
ampli.e   = etasp;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                         DSP parameters                                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Ns_bprop  = Nstep;                  % SSFM and ESSFM backpropagation steps
dsp       = DSP(ch,Ns_bprop);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                         Signal parameters                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Nsymb     = sym_length;                % number of symbols
Nt        = 2;                         % points x symbol
sig       = Signal(Nsymb,Nt,symbrate);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                         Some System Parameters                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nfft      = sig.NSYMB*sig.NT;
Loss      = 10^(-Gerbio*0.1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                         Pulse parameters                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pls.shape   = 'RRC';
pls.bw      = 1.0;                       % duty cycle
pls.ord     = 0.1;                       % pulse roll-off
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                       SNR Matched filter                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Hf_SNR        = transpose(filt(pls,sig.FN));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                     Signal Modulation (for SNR)                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[patx_tx(:,1), patmatx_tx]    = Pattern.random(4,sig.NSYMB);
[paty_tx(:,1), patmaty_tx]    = Pattern.random(4,sig.NSYMB);

E                             = Laser.GetLaserSource(1, nfft);

set(sig,'FIELDX'    ,Modulator.ApplyModulation(E, 2*patmatx_tx-1, sig, pls));
set(sig,'FIELDX_TX' ,Modulator.ApplyModulation(E, 2*patmatx_tx-1, sig, pls));
set(sig,'FIELDY'    ,Modulator.ApplyModulation(E, 2*patmaty_tx-1, sig, pls));
set(sig,'FIELDY_TX' ,Modulator.ApplyModulation(E, 2*patmaty_tx-1, sig, pls));

system.Nspan     = Nspan;
system.Loss      = Loss;
system.mfil_snr  = Hf_SNR;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

options = optimset('Display','off','TolX',3e-1);
funmax = @(Ps_dBm) -par_disp_comp_snr(Ps_dBm,ch,dsp,sig,ampli,system);

[P_max,max_snr] = fminbnd(funmax,-8,8,options);
max_snr = -max_snr;
end

