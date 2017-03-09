function [ max_snr ] = disp_comp_maxSNR(ch,dsp,sig,etasp,Nspan)
%DISP_COMP_MAXSNR Summary of this function goes here
%   Detailed explanation goes here
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                         Amplifier parameters                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Gerbio    = ch.alphalin/(log(10)*1e-4)*ch.Lf*1e-3;
ampli.G   = Gerbio;
ampli.e   = etasp;
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
%              Modulation of the Signal used for SNR  calc               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                       SNR Matched filter                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Hf_SNR        = transpose(filt(pls,sig.FN));
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

[P_max,max_snr] = fminbnd(funmax,-2,10,options);


end

