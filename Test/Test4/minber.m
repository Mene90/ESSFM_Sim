function [ min_ber ] = minber(ch,dsp,sig,t_sig,NC,etasp)
%MINBER Summary of this function goes here
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
t_nfft    = t_sig.NSYMB*t_sig.NT;
nfft      = sig.NSYMB*sig.NT;
Nspan     = 40;
Loss      = 10^(-Gerbio*0.1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                         Pulse parameters                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pls.shape   = 'RRC';
pls.bw      = 1.0;                       % duty cycle
pls.ord     = 0.2;                       % pulse roll-off
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                         Matched filter                                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Hf          = transpose(filt(pls,t_sig.FN));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                         ESSFM PARAMETERS                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
options = optimset('Algorithm','trust-region-reflective','Display','off',...
    'Jacobian','off','DerivativeCheck','off','TolFun',1e-13,'TolX',1e-13);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                         Training Signal                                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[patx(:,1), patmatx]    = Pattern.debruijn(1,4,t_sig.NSYMB);
[paty(:,1), patmaty]    = Pattern.debruijn(2,4,t_sig.NSYMB);

E                       = Laser.GetLaserSource(1, t_nfft);

set(t_sig,'FIELDX_TX',Modulator.ApplyModulation(E, 2*patmatx-1, t_sig, pls));
set(t_sig,'FIELDX'   ,Modulator.ApplyModulation(E, 2*patmatx-1, t_sig, pls));
set(t_sig,'FIELDY_TX',Modulator.ApplyModulation(E, 2*patmaty-1, t_sig, pls));
set(t_sig,'FIELDY'   ,Modulator.ApplyModulation(E, 2*patmaty-1, t_sig, pls));

essfm.NC     = NC;
essfm.opt    = options;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%              Modulation of the Signal used for BER  calc               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                       BER Matched filter                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Hf_BER        = transpose(filt(pls,sig.FN));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[patx_tx(:,1), patmatx_tx]    = Pattern.random(4,sig.NSYMB);
[paty_tx(:,1), patmaty_tx]    = Pattern.random(4,sig.NSYMB);

E                                 = Laser.GetLaserSource(1, nfft);

set(sig,'FIELDX'    ,Modulator.ApplyModulation(E, 2*patmatx_tx-1, sig, pls));
set(sig,'FIELDX_TX' ,Modulator.ApplyModulation(E, 2*patmatx_tx-1, sig, pls));
set(sig,'FIELDY'    ,Modulator.ApplyModulation(E, 2*patmaty_tx-1, sig, pls));
set(sig,'FIELDY_TX' ,Modulator.ApplyModulation(E, 2*patmaty_tx-1, sig, pls));

system.Nspan     = Nspan;
system.Loss      = Loss;
system.mfil_ber  = Hf_BER;
system.mfil      = Hf;

pat_tx.x     = patmatx_tx;
pat_tx.y     = patmaty_tx;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
options = optimset('Display','iter','TolX',1e-1);
fun2min = @(Ps_dBm) par_ber(Ps_dBm,ch,dsp,sig,t_sig,pat_tx,ampli,system,essfm);

min_ber = fminbnd(fun2min,-3,10,options);



end

