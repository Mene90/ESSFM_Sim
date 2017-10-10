function [signals,SNRdB,ch] = TestZeroDispHighPower(link,sp,signal,amp,pdbm,distribution)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                         Link parameters                                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
LL        = link.LL;              % length [m]
alphadB   = link.attenuation;     % attenuation [dB/km]
aeff      = 80;                   % effective area [um^2]
n2        = 2.5e-20;              % nonlinear index [m^2/W]
lambda    = link.lambda;          % wavelength [nm] @ dispersion
D         = 0;                    % dispersion [ps/nm/km] @ wavelength
S         = 0;                    % slope [ps/nm^2/km] @ wavelength

Ns_prop   = link.sprop;           % number of SSFM propagation step
Nspan     = link.Nspan;           % total number of amplifiers

ch        = Channel(LL,alphadB,lambda,aeff,n2,D,S,Ns_prop);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                         DSP parameters                                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Ns_bprop = sp.bprop;                % SSFM and ESSFM backpropagation steps
dsp      = DSP(ch,Ns_bprop);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                      Signal parameters                                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
symbrate = signal.symbrate;           % symbol rate             [Gbaud]
Ps_dBm   = pdbm;                      % Power vector            [dBm]
Pavg     = 10.^(0.1*(Ps_dBm-30));     % Power vector            [W]
Nsymb    = signal.nsymb;              % number of symbols
Nt       = 1;                         % points x symbol
sig      = Signal(Nsymb,Nt,symbrate,lambda,1);
set(sig,'POWER',Pavg);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                         Amplifier parameters                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Gerbio    = alphadB*LL*1e-3;
etasp     = amp.etasp;
ampli     = Ampliflat(Pavg,ch,Gerbio,etasp,'EDFA');

if(distribution == 'HG')
    [cmapx_tx] = Pattern.halfgaussian(Nsymb);
    %     [cmapy_tx] = Pattern.halfgaussian(Nsymb);
elseif(distribution == 'G')
    [cmapx_tx] = Pattern.gaussian(Nsymb);
    %     [cmapy_tx] = Pattern.gaussian(Nsymb);
end

E      = Laser.GetLaserSource(Pavg, sig,lambda,0);

set(sig,'FIELDX'    ,cmapx_tx);
set(sig,'FIELDX_TX' ,cmapx_tx);
%      set(sig,'FIELDX'    ,cmapy_tx);
%      set(sig,'FIELDY_TX' ,cmapy_tx);

gpu_propagation(ch,Nspan,ampli,sig);
dsp.backpropagation(Pavg*10^(-Gerbio*0.1),sig,Nspan,'ssfm');
signolpn = copy(sig);
dsp.nlpnmitigation(sig);

signals(1) = sig.getproperties();
signals(2) = signolpn.getproperties();
signals(1).type = 'NLPN_mitigation = true';
signals(2).type = 'NLPN_mititagion = false';

SNRdB  = 10*log10(1/symbrate/10^9/ampli.N0/Nspan);

end



