function [signals,SNRdB,ch] = DM_subcarrier_propagation( link,sp,signal,sub_signal,amp,pdbm,wdm,pls,pol,gpu )

HPLANCK = 6.62606896e-34;       % Planck's constant [J*s]
CLIGHT = 299792458;             % speed of light [m/s]

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

ch       = Channel(LL,alphadB,lambda,aeff,n2,D,S,Ns_prop);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                     Comp Link parameters                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
comp_alphadB   = 0.57;                       % attenuation [dB/km]
aeff      = 30;                              % effective area [um^2]
n2        = 4.8e-20;                         % nonlinear index [m^2/W]
lambda    = link.lambda;                     % wavelength [nm] @ dispersion
D         = -100;                            % dispersion [ps/nm/km] @ wavelength
S         = 0;                               % slope [ps/nm^2/km] @ wavelength
comp_LL   = -link.disp*link.LL/1e3/D*1e3;    % comp. fiber length [m];
Ns_prop   = 1000;                            % number of SSFM propagation step

comp_ch   = Channel(comp_LL,comp_alphadB,lambda,aeff,n2,D,S,Ns_prop);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                         DSP parameters                                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Ns_bprop      = sp.bprop;              % SSFM and ESSFM backpropagation steps
Ns_bprop_comp = 100;
dsp           = DSP(ch,Ns_bprop);
dsp_comp      = DSP(comp_ch,Ns_bprop_comp);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                      Signal parameters (wdm)                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
symbrate = signal.symbrate;           % symbol rate             [Gbaud]
Nsymb    = signal.nsymb;              % number of symbols
Nt       = signal.nt;                 % points x symbol
sig      = Signal(Nsymb,Nt,symbrate,lambda,signal.nc);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                   Signal parameters (sub carrier)                      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
symbrate = sub_signal.symbrate;           % symbol rate             [Gbaud]
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
Gerbio01  = alphadB*LL*1e-3;
Gerbio02  = comp_alphadB*comp_LL*1e-3;
etasp     = amp.etasp;
amptype   = amp.type;

Pu0       = 10.^(0.1*(pdbm-30));         % Power [W]
Pu1       = 10.^(0.1*((pdbm-4)-30));       % Power in dc [W]
Gm1       = (10^(Gerbio01*0.1)-1.d0)/Pu0 + (10^(Gerbio02*0.1)-1.d0)/Pu1 ;
ampli     = Ampliflat([],ch,Gm1,etasp,amptype,Nspan);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Subcarrier creation
set(sig,'POWER',Pu0);
sc_generation(Pu0,lambda,Nsymb,pol,pls,sub_signal,sub_sig)

%% Sub carrier Multiplexing
if (not(sub_signal.nsc == 1))
    for ii = 1:sig.NCH
        MuxDemux.Mux(sub_sig{ii}.SUB_FIELDX,sub_sig{ii}.SUB_FIELDY,sub_sig{ii});
    end
end

%% Sub carrier Modulation and channel multiplexing
Laser.GetLaserSource(Pu0,sig,lambda,0);
if (not(sub_signal.nsc == 1))
    for ii = 1:sig.NCH
        Eoptx(:,ii) = Modulator.ApplyModulation([],sub_sig{ii}.FIELDX,sig,pls);
        if(pol == 2 )
            Eopty(:,ii) = Modulator.ApplyModulation([],sub_sig{ii}.FIELDY,sig,pls);
        end
    end
else
    for ii = 1:sig.NCH
        Eoptx(:,ii) = Modulator.ApplyModulation([],sub_sig{ii}.SUB_FIELDX,sig,pls);
        if(pol == 2 )
            Eopty(:,ii) = Modulator.ApplyModulation([],sub_sig{ii}.SUB_FIELDY,sig,pls);
        end
    end
end

if(pol == 2 )
    MuxDemux.Mux(Eoptx,Eopty,sig);
else
    MuxDemux.Mux(Eoptx,[],sig);
end

%% Propagation of the multiplexed signal 
    
   set(sig,'FIELDX', gpuArray(complex(get(sig,'FIELDX'))));
   set(sig,'FIELDY', gpuArray(complex(get(sig,'FIELDY'))));

    for i = 1:Nspan
        set(sig,'POWER',Pu0);
        sing_span_propagation(ch,sig,gpu);
        set(sig,'POWER',Pu1);
        sing_span_propagation(comp_ch,sig,gpu);
    end

    AddNoise(ampli,sig);
    set(sig,'POWER', Pu0);
    set(sig,'FIELDX', gather(get(sig,'FIELDX')));
    set(sig,'FIELDY', gather(get(sig,'FIELDY')));
            
    g = gpuDevice(1);
    reset(g);

%% Demultiplexing
oHf       = myfilter(oftype,sig.FN,obw,0);    
[zfieldx,zfieldy] = MuxDemux.Demux(sig,oHf,cch);
set(sig,'FIELDX',zfieldx(:,cch));
set(sig,'FIELDY',zfieldy(:,cch));

%% Backpropagation

 %backpropagation(dsp,Pu0*10^(-Gerbio*0.1),sig,Nspan,'ssfm');
 
    set(sig,'FIELDX', gpuArray(complex(get(sig,'FIELDX'))));
    set(sig,'FIELDY', gpuArray(complex(get(sig,'FIELDY'))));
 
    for i = 1:Nspan
        set(sig,'POWER',Pu1);
        backpropagation(dsp_comp,Pu1*10^(-Gerbio02*0.1),sig,1,'ssfm',0);
        set(sig,'POWER',Pu0);
        backpropagation(dsp,Pu0*10^(-Gerbio01*0.1),sig,1,'ssfm',0);        
     end
    
    set(sig,'FIELDX', gather(get(sig,'FIELDX')));
    set(sig,'FIELDY', gather(get(sig,'FIELDY')));
    
    g      = gpuDevice(1);
    reset(g);

%% Subcarrier Demultiplexing
 
 dsp.downsampling(sig);
 set(sig, 'FIELDX_TX', sub_sig{cch}.FIELDX_TX); 
 set(sig, 'FIELDY_TX', sub_sig{cch}.FIELDY_TX);
 
 if (not(sub_signal.nsc == 1))
     oHf      = myfilter(oftype,sub_sig{cch}.FN,obw,0);
     set(sig,'FN',sub_sig{cch}.FN);
     set(sig,'SYMBOLRATE',sub_signal.symbrate);
     set(sig,'NCH',sub_signal.nsc)
     set(sig, 'FIELDX_TX', sub_sig{cch}.FIELDX_TX);
     set(sig, 'FIELDY_TX', sub_sig{cch}.FIELDY_TX);
     set(sig, 'NT', sub_signal.nt);
     Laser.GetLaserSource(Pu0,sig,lambda,0); 
     if(pol == 2)
         [sig.SUB_FIELDX,sig.SUB_FIELDY] = MuxDemux.Demux(sig,oHf,0);
     else
         [sig.SUB_FIELDX] = MuxDemux.Demux(sig,oHf,0);
     end
     dsp.scdownsampling(sig);
 end
 
 signals      = sig.getproperties();
 SNRdB        = 10*log10(1/(symbrate*sub_signal.nsc)/10^9/amp.N0);
 
end

function sc_generation(Pu0,lambda,Nsymb,pol,pls,sub_signal,sub_sig)

for ii = 1:sig.NCH
    set(sub_sig{ii},'POWER',Pu0);
    Laser.GetLaserSource(Pu0,sub_sig{ii},lambda,0);
    for j = 1:sub_signal.nsc
        [sub_cmapx(:,j)] = Pattern.subc_gaussian(Nsymb,1/sqrt(sub_signal.nsc));
        if(pol == 2)
            [sub_cmapy(:,j)] = Pattern.subc_gaussian(Nsymb,1/sqrt(sub_signal.nsc));
        end
        if (not(sub_signal.nsc == 1))
            sub_Eoptx(:,j)   = Modulator.ApplyModulation([],sub_cmapx(:,j),sub_sig{ii},pls);
            if(pol == 2 )
                sub_Eopty(:,j) = Modulator.ApplyModulation([],sub_cmapy(:,j),sub_sig{ii},pls);
            end
        else
            sub_Eoptx(:,j)   = sub_cmapx(:,j);
            if (pol == 2 )
                sub_Eopty(:,j) = sub_cmapy(:,j);
            end
        end
    end
    
    set(sub_sig{ii},'FIELDX_TX' ,sub_cmapx);
    set(sub_sig{ii},'SUB_FIELDX',sub_Eoptx); 
    if (pol == 2)
        set(sub_sig{ii},'FIELDY_TX' ,sub_cmapy);
        set(sub_sig{ii},'SUB_FIELDY',sub_Eopty);
    end
    
%     sub_carriers{ii} = sub_sig{ii}.getproperties();
end
end

