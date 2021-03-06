function [signals,SNRdB,ch] = Test_mux(link,sp,signal,amp,pdbm,wdm,pls,gpu)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                         Link parameters                                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
LL        = link.LL;              % length [m]
alphadB   = link.attenuation;     % attenuation [dB/km]
aeff      = 80;                   % effective area [um^2]
n2        = link.nlindex;         % nonlinear index [m^2/W]
lambda    = link.lambda;          % wavelength [nm] @ dispersion
D         = link.disp;            % dispersion [ps/nm/km] @ wavelength
S         = 0;                    % slope [ps/nm^2/km] @ wavelength

Ns_prop   = link.sprop;           % number of SSFM propagation step
Nspan     = link.Nspan;           % total number of amplifiers

pmd       = false;                % pmd enable/disable

ch       = Channel(LL,alphadB,lambda,aeff,n2,D,S,Ns_prop);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                         DSP parameters                                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Ns_bprop = sp.bprop;                  % SSFM and ESSFM backpropagation steps
dsp      = DSP(ch,Ns_bprop);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                      Signal parameters                                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
symbrate = signal.symbrate;           % symbol rate             [Gbaud]
Ps_dBm   = pdbm;                      % Power vector            [dBm]
Pavg     = 10.^(0.1*(Ps_dBm -30));    % Power vector            [W]
Plen     = length(Ps_dBm);  
Nsymb    = signal.nsymb;              % number of symbols
Nt       = signal.nt;                 % points x symbol
sig      = Signal(Nsymb,Nt,symbrate,lambda,signal.nc);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                         WDM parameters                                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
cch     =  wdm.cch;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                         Pulse parameters                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% pls.shape   = 'RRC';                      % Shape type
% pls.bw      = 1.0;                        % duty cycle
% pls.ord     = 0.1;                        % pulse roll-off
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                         optical filter parameters                      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
oftype = 'ideal'; % optical filter type
obw    = 1;       % optical filter bandwidth 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                         Amplifier parameters                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Gerbio    = alphadB*LL*1e-3;
etasp     = amp.etasp;
amptype   = amp.type;
% ampli     = Ampliflat(Pavg,ch,Gerbio,etasp);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                         Matched filter                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Hf        = filt(pls,sig.FN);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                         Optical filter                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
oHf       = myfilter(oftype,sig.FN,obw,0);      % Remember that the in the lowpass
                                                % equivalent domain, the 3 dB bandwidth
                                                % goes from -bw/2 to + bw/2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% for i = 1:Plen
       
    set(sig,'POWER',Pavg);
    
    ampli  = Ampliflat(Pavg,ch,Gerbio,etasp,amptype,Nspan);
    E      = Laser.GetLaserSource(Pavg,sig,lambda,0); %0.400835
    
    for ii = 1:sig.NCH
        
        [cmapx(:,ii)]   = Pattern.gaussian(Nsymb);
        Eoptx(:,ii)     = Modulator.ApplyModulation(E,cmapx(:,ii),sig,pls);
        
    end
    
    set(sig,'FIELDX_TX',cmapx);
    MuxDemux.Mux(Eoptx,[],sig);
    
   
    if gpu
        set(sig,'FIELDX', gpuArray(complex(get(sig,'FIELDX'))));
        set(sig,'FIELDY', gpuArray(complex(get(sig,'FIELDY'))));
    end
    
    for i = 1:Nspan
        sing_span_propagation(ch,sig,gpu);
        AddNoise(ampli,sig);
    end
    
     if gpu
        set(sig,'FIELDX', gather(get(sig,'FIELDX')));
        set(sig,'FIELDY', gather(get(sig,'FIELDY')));
     end
    

    [zfieldx] = MuxDemux.Demux(sig,oHf,0);
    
    set(sig,'FIELDX',zfieldx(:,cch));
    set(sig,'FIELDX_TX',sig.FIELDX_TX(:,cch));
    
    if (strcmp(amptype,'Raman'))
        dsp.backpropagation(Pavg,sig,Nspan,'ssfm',gpu);
    else
        dsp.backpropagation(Pavg*10^(-Gerbio*0.1),sig,Nspan,'ssfm',gpu);
    end

%    dsp.matchedfilter(sig,Hf);
    dsp.downsampling(sig);
%    dsp.nlpnmitigation(sig);
%     
%    avgber       = ErrorEstimation.BER(sig);
    signals      = sig.getproperties();
    SNRdB        = 10*log10(1/symbrate/10^9/ampli.N0);
    
end