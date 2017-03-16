function [ max_snr,P_max ] = ESSFM_MAX_SNR( Nstep,NC,sym_length,n_prop_steps,etasp,R,Nspan)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                      Global Signal parameters                          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
symbrate  = R;                    % symbol rate [Gbaud]
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
%                      Trainin Signal parameters                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Nsymb     = 2^10;                % number of symbols
Nt        = 2;                   % points x symbol
t_sig     = Signal(Nsymb,Nt,symbrate);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                            Signal parameters                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Nsymb     = sym_length;                % number of symbols
Nt        = 2;                         % points x symbol
sig       = Signal(Nsymb,Nt,symbrate);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           System Parameters                            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t_nfft    = t_sig.NSYMB*t_sig.NT;
nfft      = sig.NSYMB*sig.NT;
Loss      = 10^(-Gerbio*0.1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                         ESSFM PARAMETERS                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
options = optimset('Algorithm','trust-region-reflective','Display','off',...
    'Jacobian','off','DerivativeCheck','off','TolFun',1e-13,'TolX',1e-13);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                         Pulse parameters                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pls.shape   = 'RRC';
pls.bw      = 1.0;                       % duty cycle
pls.ord     = 0.1;                       % pulse roll-off
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                         Matched filter                                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Hf          = transpose(filt(pls,t_sig.FN));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                         Training Signal Modulation                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
seedx = 1;
seedy = t_sig.NSYMB/2^3;
[patx(:,1), patmatx]    = Pattern.debruijn(seedx,4,t_sig.NSYMB);
[paty(:,1), patmaty]    = Pattern.debruijn(seedy,4,t_sig.NSYMB);

E                       = Laser.GetLaserSource(1, t_nfft);

set(t_sig,'FIELDX_TX',1./sqrt(2.)*((2*patmatx(:,1)-1)+1i*(2.*patmatx(:,2)-1)));
set(t_sig,'FIELDX'   ,Modulator.ApplyModulation(E, 2*patmatx-1, t_sig, pls));
set(t_sig,'FIELDY_TX',1./sqrt(2.)*((2*patmaty(:,1)-1)+1i*(2.*patmaty(:,2)-1)));
set(t_sig,'FIELDY'   ,Modulator.ApplyModulation(E, 2*patmaty-1, t_sig, pls));

essfm.NC     = NC;
essfm.opt    = options;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                       SNR Matched filter                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Hf_SNR        = transpose(filt(pls,sig.FN));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                         Signal Modulation (for SNR)                    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[patx_tx(:,1), patmatx_tx]    = Pattern.random(4,sig.NSYMB);
[paty_tx(:,1), patmaty_tx]    = Pattern.random(4,sig.NSYMB);

E                             = Laser.GetLaserSource(1, nfft);

set(sig,'FIELDX'    ,Modulator.ApplyModulation(E, 2*patmatx_tx-1, sig, pls));
set(sig,'FIELDX_TX' ,1./sqrt(2.)*((2*patmatx_tx(:,1)-1)+1i*(2.*patmatx_tx(:,2)-1)));
set(sig,'FIELDY'    ,Modulator.ApplyModulation(E, 2*patmaty_tx-1, sig, pls));
set(sig,'FIELDY_TX' ,1./sqrt(2.)*((2*patmaty_tx(:,1)-1)+1i*(2.*patmaty_tx(:,2)-1)));

system.Nspan     = Nspan;
system.Loss      = Loss;
system.mfil_snr  = Hf_SNR;
system.mfil      = Hf;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
options = optimset('Display','off','TolX',1e-2);
funmax = @(Ps_dBm) -par_essfm_snr(Ps_dBm,ch,dsp,sig,t_sig,ampli,system,essfm);

[P_max,max_snr] = fminbnd(funmax,-2,10,options);
max_snr = -max_snr;
end

