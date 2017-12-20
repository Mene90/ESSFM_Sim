function [ signals,SNRdB,ch ] = TestDispOrKerrComp( link,sp,signal,amp,pdbm,distribution,compensation,gpu,pls)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                         Link parameters                                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
LL        = link.LL;              % length [m]
alphadB   = link.attenuation;     % attenuation [dB/km]
aeff      = 80;                   % effective area [um^2]
n2        = link.n2;              % nonlinear index [m^2/W]
lambda    = link.lambda;          % wavelength [nm] @ dispersion
D         = link.disp;            % dispersion [ps/nm/km] @ wavelength
S         = 0;                    % slope [ps/nm^2/km] @ wavelength

Ns_prop   = link.sprop;           % number of SSFM propagation step
Nspan     = link.Nspan;           % total number of amplifiers

ch        = Channel(LL,alphadB,lambda,aeff,n2,D,S,Ns_prop);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                     Comp Link parameters                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
comp_alphadB   = 0;                     % attenuation [dB/km]
aeff      = 30;                         % effective area [um^2]
n2        = 0;                          % nonlinear index [m^2/W]
lambda    = link.lambda;                % wavelength [nm] @ dispersion
D         = -100;                       % dispersion [ps/nm/km] @ wavelength
S         = 0;                          % slope [ps/nm^2/km] @ wavelength
comp_LL   = -link.disp*link.LL/1e3/D*1e3;    % comp. fiber length [m];           

Ns_prop   = 1;                           % number of SSFM propagation step

comp_ch   = Channel(comp_LL,comp_alphadB,lambda,aeff,n2,D,S,Ns_prop);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                         DSP parameters                                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Ns_bprop = sp.bprop;                % SSFM and ESSFM backpropagation steps
Ns_bprop_comp = 1;
dsp       = DSP(ch,Ns_bprop);
dsp_dbp   = DSP(ch,link.sprop);
dsp_comp  = DSP(comp_ch,Ns_bprop_comp);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                      Signal parameters                                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
symbrate = signal.symbrate;           % symbol rate             [Gbaud]
Ps_dBm   = pdbm;                      % Power vector            [dBm]
Pavg     = 10.^(0.1*(Ps_dBm-30));     % Power vector            [W]
Nsymb    = signal.nsymb;              % number of symbols
Nt       = signal.nt;                 % points x symbol
sig      = Signal(Nsymb,Nt,symbrate,lambda,1);
set(sig,'POWER',Pavg);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                         Amplifier parameters                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Gerbio    = alphadB*LL*1e-3;
Gerbio    = alphadB*LL*1e-3+comp_alphadB*comp_LL*1e-3;
etasp     = amp.etasp;
amptype   = amp.type;
ampli     = Ampliflat(Pavg,ch,Gerbio,etasp,amptype,Nspan);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                         Matched filter                                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Hf        = filt(pls,sig.FN);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

E      = Laser.GetLaserSource(Pavg,sig,lambda);

if(distribution == 'HG')
    [cmapx_tx]  = Pattern.halfgaussian(Nsymb);
    Eoptx       = Modulator.ApplyModulation(E,cmapx_tx,sig,pls);
elseif(distribution == 'G')
    [cmapx_tx] = Pattern.gaussian(Nsymb);
    Eoptx      = Modulator.ApplyModulation(E,cmapx_tx,sig,pls);
end

% set(sig,'FIELDX'    ,Eoptx);
% set(sig,'FIELDX_TX' ,Eoptx);
set(sig,'FIELDX'    ,Eoptx);
set(sig,'FIELDX_TX' ,cmapx_tx);

if (strcmp(compensation,'inline'))
    
    if gpu
        set(sig,'FIELDX', gpuArray(complex(get(sig,'FIELDX'))));
        set(sig,'FIELDY', gpuArray(complex(get(sig,'FIELDY'))));
    end
    
    %propagation(ch,Nspan,ampli,sig);
    
    for i = 1:Nspan
        sing_span_propagation(ch,sig,gpu);
        sing_span_propagation(comp_ch,sig,gpu);
        AddNoise(ampli,sig);
    end
    
    if gpu
        set(sig,'FIELDX', gather(get(sig,'FIELDX')));
        set(sig,'FIELDY', gather(get(sig,'FIELDY')));
    end
    
    %backpropagation(dsp,Pavg*10^(-Gerbio*0.1),sig,Nspan,'ssfm');
    
    sig_dbp = copy(sig);
    
    for i = 1:Nspan
        backpropagation(dsp_comp,Pavg*10^(-Gerbio*0.1),sig_dbp,1,'ssfm',gpu);
        backpropagation(dsp_dbp,Pavg*10^(-Gerbio*0.1),sig_dbp,1,'ssfm',gpu);
    end
   
else 
    
    if gpu
        gpu_propagation(ch,Nspan,ampli,sig);
    else
        propagation(ch,Nspan,ampli,sig);
    end
    
    sig_dbp = copy(sig);
    
    backpropagation(dsp,Pavg*10^(-Gerbio*0.1),sig,Nspan,compensation,gpu);
    backpropagation(dsp_dbp,Pavg*10^(-Gerbio*0.1),sig_dbp,Nspan,'ssfm',gpu);

end

dsp.matchedfilter(sig,Hf);
dsp.downsampling(sig);

dsp.matchedfilter(sig_dbp,Hf);
dsp.downsampling(sig_dbp);

signals     = sig.getproperties();
signals_dbp = sig_dbp.getproperties();

SNRdB  = 10*log10(1/symbrate/10^9/ampli.N0/Nspan);


end

