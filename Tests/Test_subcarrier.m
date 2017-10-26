function [signals,SNRdB,ch] = Test_subcarrier(link,sp,signal,sub_signal,amp,pdbm,wdm,pls,Nsc)
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
%                      Signal parameters (wdm)                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
symbrate = signal.symbrate;           % symbol rate             [Gbaud]
Ps_dBm   = pdbm;                      % Power vector            [dBm]
Pavg     = 10.^(0.1*(Ps_dBm-30));    % Power vector            [W]
Plen     = length(Ps_dBm);  
Nsymb    = signal.nsymb;              % number of symbols
Nt       = signal.nt;                 % points x symbol
sig      = Signal(Nsymb,Nt,symbrate,lambda,signal.nc);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                   Signal parameters (sub carrier)                      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
symbrate = sub_signal.symbrate;       % symbol rate             [Gbaud]
Ps_dBm   = pdbm;                      % Power vector            [dBm]
Pavg     = 10.^(0.1*(Ps_dBm -30));    % Power vector            [W]
Plen     = length(Ps_dBm);  
Nsymb    = sub_signal.nsymb;              % number of symbols
Nt       = sub_signal.nt;                 % points x symbol
for ii = 1:signal.nc
    sub_sig{ii}  = Signal(Nsymb,Nt,symbrate,lambda,sub_signal.nsc);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                         WDM parameters                                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
cch     =  wdm.cch;
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

set(sig,'POWER',Pavg);


ampli  = Ampliflat(Pavg,ch,Gerbio,etasp,amptype,Nspan);

%% Subcarrier creation
for ii = 1:sig.NCH
    set(sub_sig{ii},'POWER',Pavg);
    Laser.GetLaserSource(Pavg,sub_sig{ii},lambda,0.10017);
    
    for j = 1:Nsc
        [sub_cmapx(:,j)]  = Pattern.subc_gaussian(Nsymb,1/sqrt(sub_signal.nsc));
         sub_Eoptx(:,j)   = Modulator.ApplyModulation([],sub_cmapx(:,j),sub_sig{ii},pls);
    end
%     sc_Eoptx{ii}    = sub_Eoptx;
    set(sub_sig{ii},'FIELDX_TX' ,sub_cmapx);
    set(sub_sig{ii},'SUB_FIELDX',sub_Eoptx);   
end

%% Sub carrier Multiplexing
for ii = 1:sig.NCH
    MuxDemux.Mux(sub_sig{ii}.SUB_FIELDX,[],sub_sig{ii},2);
end

%% Channel Multiplexing
Laser.GetLaserSource(Pavg,sig,lambda,0.400835);
for ii = 1:sig.NCH    
    Eoptx(:,ii) = Modulator.ApplyModulation([],sub_sig{ii}.FIELDX,sig,pls);
end
MuxDemux.Mux(Eoptx,[],sig,0);

%% Propagation
    for i = 1:Nspan
        sing_span_propagation(ch,sig,'true')
    end
    AddNoise(ampli,sig);

%% Channel Demultiplexing
oHf       = myfilter(oftype,sig.FN,obw,0);    
[zfieldx] = MuxDemux.Demux(sig,oHf,cch,0);
set(sig,'FIELDX',zfieldx(:,cch));

%% Backpropagation
% 
    if (strcmp(amptype,'Raman'))
        dsp.backpropagation(Pavg,sig,Nspan,'ssfm',1);
    else
        dsp.backpropagation(Pavg*10^(-Gerbio*0.1),sig,Nspan,'ssfm',1);
    end

%% Subcarrier Demultiplexing
 dsp.downsampling(sig); 
 oHf      = myfilter(oftype,sub_sig{1}.FN,obw,0);
 set(sig,'FN',sub_sig{1}.FN);
 set(sig,'SYMBOLRATE',sub_signal.symbrate); 
 set(sig,'NCH',sub_signal.nsc)
 set(sig, 'FIELDX_TX', sub_sig{cch}.FIELDX_TX);
 set(sig, 'NT', sub_signal.nt);
 Laser.GetLaserSource(Pavg,sig,lambda,0.10017);
 [sig.SUB_FIELDX] = MuxDemux.Demux(sig,oHf,0,2);
 dsp.scdownsampling(sig);
 signals      = sig.getproperties();
 SNRdB        = 10*log10(1/symbrate/10^9/ampli.N0);
 
end

