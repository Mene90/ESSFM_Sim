%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                         Link parameters                                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
LL        = 1e5;                  % length [m]
alphadB   = 0.2;                  % attenuation [dB/km]
aeff      = 80;                   % effective area [um^2]
n2        = 2.5e-20;              % nonlinear index [m^2/W]
lambda    = 1550;                 % wavelength [nm] @ dispersion
D         = 17;                   % dispersion [ps/nm/km] @ wavelength
S         = 0;                    % slope [ps/nm^2/km] @ wavelength

Ns_prop   = 10;                   % number of SSFM propagation step
Nspan     = 40;                   % total number of amplifiers

pmd       = false;                % pmd enable/disable

ch        = Channel(LL,alphadB,lambda,aeff,n2,D,S,Ns_prop);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                         DSP parameters                                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Ns_bprop = 1;                  % SSFM and ESSFM backpropagation steps
dsp      = DSP(ch,Ns_bprop);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                      Signal parameters                                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
symbrate = 32;                        % symbol rate             [Gbaud]
Ps_dBm   = -1;                         % Power vector            [dBm]
Pavg     = 10.^(0.1*(Ps_dBm-30));     % Power vector            [W]
Plen     = length(Ps_dBm);  
Nsymb    = 2^16;                      % number of symbols
Nt       = 2;                         % points x symbol
sig      = Signal(Nsymb,Nt,symbrate);
set(sig,'POWER',Pavg);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                         Pulse parameters                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pls.shape   = 'RRC';                     % Shape type
pls.bw      = 1.0;                       % duty cycle
pls.ord     = 0.2;                       % pulse roll-off
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                         Amplifier parameters                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Gerbio    = alphadB*LL*1e-3;
etasp     = 2;
ampli     = Ampliflat(Pavg,ch,Gerbio,etasp);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                         Matched filter                                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Hf        = transpose(filt(pls,sig.FN));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[patx_tx(:,1), patmatx_tx] = Pattern.random(4,Nsymb);
[paty_tx(:,1), patmaty_tx] = Pattern.random(4,Nsymb);
E                          = Laser.GetLaserSource(Pavg, Nsymb*Nt);

set(sig,'FIELDX'    ,Modulator.ApplyModulation(E, 2*patmatx_tx-1, sig, pls));
set(sig,'FIELDY'    ,Modulator.ApplyModulation(E, 2*patmaty_tx-1, sig, pls));
set(sig,'FIELDX_TX' ,Modulator.ApplyModulation(E, 2*patmatx_tx-1, sig, pls));
set(sig,'FIELDY_TX' ,Modulator.ApplyModulation(E, 2*patmaty_tx-1, sig, pls));

ch.propagate(Nspan,ampli,sig);

dsp.setessfmcoeff(5,Pavg,2^10,2,symbrate,pls,ampli,Nspan);
dsp.backpropagation(Pavg*10^(-alphadB*LL*1e-3*0.1),sig,Nspan,'essfm');
dsp.matchedfilter(sig,Hf);
dsp.nlpnmitigation(sig);

avgber = Ber.BerEstimation(sig);

