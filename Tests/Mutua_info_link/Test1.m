function [signals,SNRdB,ch] = Test1(nlindex,disp,link,sp,signal,amp,pdbm)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                         Link parameters                                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
LL        = link.LL;              % length [m]
alphadB   = link.attenuation;     % attenuation [dB/km]
aeff      = 80;                   % effective area [um^2]
n2        = nlindex;%2.5e-20;     % nonlinear index [m^2/W]
lambda    = link.lambda;          % wavelength [nm] @ dispersion
D         = disp;%17;             % dispersion [ps/nm/km] @ wavelength
S         = 0;                    % slope [ps/nm^2/km] @ wavelength

Ns_prop   = link.sprop;           % number of SSFM propagation step
Nspan     = link.Nspan;           % total number of amplifiers

ch        = Channel(LL,alphadB,lambda,aeff,n2,D,S,Ns_prop);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                     Comp Link parameters                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% comp_alphadB   = 9e-4;                  % attenuation [dB/km]
% aeff      = 20;                         % effective area [um^2]
% n2        = 2.5e-20;                    % nonlinear index [m^2/W]
% lambda    = link.lambda;                % wavelength [nm] @ dispersion
% D         = -1010;                      % dispersion [ps/nm/km] @ wavelength
% S         = 0;                          % slope [ps/nm^2/km] @ wavelength
% comp_LL        = -disp*link.LL/1e3/D*1e3;    % comp. fiber length [m];           
% 
% Ns_prop   = link.sprop/10;                 % number of SSFM propagation step
% 
% comp_ch   = Channel(comp_LL,comp_alphadB,lambda,aeff,n2,D,S,Ns_prop);
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
% Plen     = length(Ps_dBm);  
Nsymb    = signal.nsymb;              % number of symbols
Nt       = 1;                         % points x symbol
sig      = Signal(Nsymb,Nt,symbrate,lambda,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                         Pulse parameters                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pls.shape   = 'RRC';                     % Shape type
pls.bw      = 1.0;                       % duty cycle
pls.ord     = 0.2;                       % pulse roll-off
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                         Amplifier parameters                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Loss        = alphadB*LL*1e-3;
Gerbio    = alphadB*LL*1e-3;
% Gerbio    = alphadB*LL*1e-3+comp_alphadB*comp_LL*1e-3;
etasp     = amp.etasp;
% ampli     = Ampliflat(Pavg,ch,Gerbio,etasp);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                         Matched filter                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Hf        = filt(pls,sig.FN);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% for i = 1:Plen
       
    set(sig,'POWER',Pavg);
    
    ampli  = Ampliflat(Pavg,ch,Gerbio,etasp);   
    
    [cmapx_tx] = Pattern.halfgaussian(Nsymb);
%     [cmapy_tx] = Pattern.gaussian(Nsymb);
    E      = Laser.GetLaserSource(Pavg, sig,lambda,0);
%     pat   = pat2stars(patx,'qpsk')*exp(1i*pi/4);
%     patmat_tx = 2*patmatx_tx-1;
%     patx   = 1./sqrt(2.)*(patmat_tx(:,1)+1i*patmat_tx(:,2));
     set(sig,'FIELDX'    ,cmapx_tx);
%      set(sig,'FIELDX'    ,Modulator.ApplyModulation(E,cmapx_tx,sig,pls));
%      set(sig,'FIELDY'    ,Modulator.ApplyModulation(E,cmapy_tx,sig,pls));
     set(sig,'FIELDX_TX' ,cmapx_tx);
%      set(sig,'FIELDY_TX' ,cmapy_tx);
%      
%      set(sig,'FIELDX', gpuArray(complex(get(sig,'FIELDX'))));
%      set(sig,'FIELDY', gpuArray(complex(get(sig,'FIELDY'))));
     
     propagate(ch,Nspan,ampli,sig);
%     
%     for j=1:Nspan
%         AddNoise(ampli,sig);
%         ssfm(ch,get(sig,'POWER'),sig);
%         ssfm(comp_ch,get(sig,'POWER')*10^(-Loss*0.1),sig);
%     end
% %                 
%     set(sig,'FIELDX', gather(get(sig,'FIELDX')));
%     set(sig,'FIELDY', gather(get(sig,'FIELDY')));
%     
    dsp.backpropagation(Pavg*10^(-Gerbio*0.1),sig,Nspan,'ssfm');
%     dsp.matchedfilter(sig,Hf);
    signolpn = copy(sig);
    dsp.nlpnmitigation(sig);
    
%      set(sig,'FIELDX',sig.FIELDX(1:sig.NT:end));
%      set(sig,'FIELDY',sig.FIELDY(1:sig.NT:end));
%     avgber    = ErrorEstimation.BER(sig);
%     mse        = ErrorEstimation.MSE(sig);
    signals(1) = sig.getproperties();   
    signals(2) = signolpn.getproperties();
    signals(1).type = 'NLPN_mitigation = true';
    signals(2).type = 'NLPN_mititagion = false';
    SNRdB  = 10*log10(1/symbrate/10^9/ampli.N0/Nspan);
%     sig          = Signal(Nsymb,Nt,symbrate,lambda,1); 
    
% end
%         gamma   = ch.gamma;
%         b2      = ch.b2;
%         alphadB = ch.alphadB;
%         savefile =strcat('TestResult_gamma_',num2str(round(ch.gamma*10^3)),'_beta_',num2str(round(ch.b2*10^28)));
%         save(savefile,'signals','gamma','b2','alphadB','SNRdB');
end

