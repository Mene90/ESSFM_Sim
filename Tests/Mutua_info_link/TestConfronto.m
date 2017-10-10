function [signals,SNRdB,ch] = TestConfronto(link,sp,signal,amp,pdbm,singlepol,gpu)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                         Link parameters                                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
LL        = link.LL;              % length [m]
alphadB   = link.attenuation;     % attenuation [dB/km]
aeff      = 80;                   % effective area [um^2]
n2        = link.nlindex;%2.5e-20;% nonlinear index [m^2/W]
lambda    = link.lambda;          % wavelength [nm] @ dispersion
D         = link.disp;%17;        % dispersion [ps/nm/km] @ wavelength
S         = 0;                    % slope [ps/nm^2/km] @ wavelength

Ns_prop   = link.sprop;           % number of SSFM propagation step
Nspan     = link.Nspan;           % total number of amplifiers

ch        = Channel(LL,alphadB,lambda,aeff,n2,D,S,Ns_prop);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                     Comp Link parameters                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
comp_alphadB   = 0;                     % attenuation [dB/km]
aeff      = 20;                         % effective area [um^2]
n2        = 0;                          % nonlinear index [m^2/W]
lambda    = link.lambda;                % wavelength [nm] @ dispersion
D         = -100;                       % dispersion [ps/nm/km] @ wavelength
S         = 0;                          % slope [ps/nm^2/km] @ wavelength
comp_LL   = -link.disp*link.LL/1e3/D*1e3;    % comp. fiber length [m];           

Ns_prop   = 4;                           % number of SSFM propagation step

comp_ch   = Channel(comp_LL,comp_alphadB,lambda,aeff,n2,D,S,Ns_prop);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                         DSP parameters                                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Ns_bprop = sp.bprop;                % SSFM and ESSFM backpropagation steps
Ns_bprop_comp = 4;
dsp       = DSP(ch,Ns_bprop);
dsp_comp  = DSP(comp_ch,Ns_bprop_comp);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                      Signal parameters                                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
symbrate = signal.symbrate;           % symbol rate             [Gbaud]
Ps_dBm   = pdbm;                      % Power vector            [dBm]
Pavg     = 10.^(0.1*(Ps_dBm-30));     % Power vector            [W]
Nsymb    = signal.nsymb;              % number of symbols
Nt       = 2;                         % points x symbol
sig      = Signal(Nsymb,Nt,symbrate,lambda,1);
set(sig,'POWER',Pavg);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                         Pulse parameters                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pls.shape   = 'RRC';                     % Shape type
pls.bw      = 1.0;                       % duty cycle
pls.ord     = 0.25;                      % pulse roll-off
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                         Amplifier parameters                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Loss    = alphadB*LL*1e-3;
% Gerbio  = alphadB*LL*1e-3;
Gerbio    = alphadB*LL*1e-3+comp_alphadB*comp_LL*1e-3;
etasp     = amp.etasp;
ampli     = Ampliflat(Pavg,ch,Gerbio,etasp,'EDFA');   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                         Matched filter                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Hf        = filt(pls,sig.FN);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    E      = Laser.GetLaserSource(Pavg, sig,lambda,0);
    
    [cmapx_tx] = Pattern.qam64(Nsymb);
    set(sig,'FIELDX'    ,Modulator.ApplyModulation(E,cmapx_tx,sig,pls));
    set(sig,'FIELDX_TX' ,cmapx_tx);
    
    if(not(singlepol))
        [cmapy_tx] = Pattern.gaussian(Nsymb);
        set(sig,'FIELDY'    ,Modulator.ApplyModulation(E,cmapy_tx,sig,pls));
        set(sig,'FIELDY_TX' ,cmapy_tx);
    end
    
%     set(sig,'FIELDX', gpuArray(complex(get(sig,'FIELDX'))));
%     set(sig,'FIELDY', gpuArray(complex(get(sig,'FIELDY'))));
    
    %propagation(ch,Nspan,ampli,sig);
    for i = 1:Nspan
        AddNoise(ampli,sig);
        sing_span_propagation(ch,sig,gpu);
        sing_span_propagation(comp_ch,sig,gpu);
    end
        
    %backpropagation(dsp,Pavg*10^(-Gerbio*0.1),sig,Nspan,'ssfm');
    for i = 1:Nspan
        backpropagation(dsp_comp,Pavg*10^(-Gerbio*0.1),sig,1,'ssfm',gpu);
        backpropagation(dsp,Pavg*10^(-Gerbio*0.1),sig,1,'ssfm',gpu);        
    end
    
%     set(sig,'FIELDX', gather(get(sig,'FIELDX')));
%     set(sig,'FIELDY', gather(get(sig,'FIELDY')));
    
    matchedfilter(dsp,sig,Hf);
    nlpnmitigation(dsp,sig);
    
    signals = sig.getproperties();   

    SNRdB  = 10*log10(1/symbrate/10^9/ampli.N0/Nspan);
    
end

