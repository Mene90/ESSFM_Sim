function [signals,SNRdB,ch] = Test_mux(nlindex,disp,link,sp,signal,amp)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                         Link parameters                                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
LL        = link.LL;              % length [m]
alphadB   = link.attenuation;     % attenuation [dB/km]
aeff      = 80;                   % effective area [um^2]
n2        = nlindex;     % nonlinear index [m^2/W]
lambda    = link.lambda;          % wavelength [nm] @ dispersion
D         = disp;             % dispersion [ps/nm/km] @ wavelength
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
Ps_dBm   = signal.pdbm;               % Power vector            [dBm]
Pavg     = 10.^(0.1*(Ps_dBm-30));     % Power vector            [W]
Plen     = length(Ps_dBm);  
Nsymb    = signal.nsymb;              % number of symbols
Nt       = 10;                        % points x symbol
sig      = Signal(Nsymb,Nt,symbrate,lambda,5);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                         Pulse parameters                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pls.shape   = 'RRC';                     % Shape type
pls.bw      = 1.0;                       % duty cycle
pls.ord     = 0.2;                       % pulse roll-off
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                         optical filter parameters                      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
oftype = 'ideal'; % optical filter type
obw    = 2.5;     % optical filter bandwidth 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                         Amplifier parameters                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Gerbio    = alphadB*LL*1e-3;
etasp     = 2;
% ampli     = Ampliflat(Pavg,ch,Gerbio,etasp);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                         Matched filter                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Hf        = filt(pls,sig.FN);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                         Optical filter                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
oHf       = myfilter(oftype,sig.FN,obw*0.5,0);  % Remember that the in the lowpass
                                                % equivalent domain, the 3 dB bandwidth
                                                % goes from -bw/2 to + bw/2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 1:Plen
       
    set(sig,'POWER',Pavg(i));
    
    ampli  = Ampliflat(Pavg(i),ch,Gerbio,etasp);
    E     = Laser.GetLaserSource(Pavg(i),sig,lambda,0.8);
    
    for ii = 1:sig.NCH
        [cmapx(:,ii)] = Pattern.gaussian(Nsymb);        
        Eoptx(:,ii) = Modulator.ApplyModulation(E,cmapx(:,ii),sig,pls);
    end
    set(sig,'FIELDX_TX',cmapx);
    MuxDemux.Mux(Eoptx,[],sig);
    
    ch.propagate(Nspan,ampli,sig);
    [zfieldx] = MuxDemux.Demux(sig,oHf);
    set(sig,'FIELDX',zfieldx(:,1));
    
    dsp.backpropagation(Pavg(i)*10^(-Gerbio*0.1),sig,Nspan,'ssfm');
    dsp.matchedfilter(sig,Hf);
    set(sig,'FIELDX_TX',sig.FIELDX_TX(:,1));
    dsp.nlpnmitigation(sig);
%      set(sig,'FIELDX',sig.FIELDX(1:sig.NT:end));
%      set(sig,'FIELDY',sig.FIELDY(1:sig.NT:end));
    avgber       = Ber.BerEstimation(sig);
    signals{i}   = sig.getproperties();
    SNRdB{i}     = 10*log10(1/symbrate/10^9/ampli.N0/Nspan);
    sig          = Signal(Nsymb,Nt,symbrate,lambda,5); 
    
end
end