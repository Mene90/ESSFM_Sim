function [signals,SNRdB,ch] = Test_subcarrier_CompLink( link,sp,signal,sub_signal,amp,pdbm,wdm,pls,pol,gpu )

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

pmd       = false;                % pmd enable/disable

ch       = Channel(LL,alphadB,lambda,aeff,n2,D,S,Ns_prop);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                     Comp Link parameters                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
comp_alphadB   = 0.57;                  % attenuation [dB/km]
aeff      = 30;                         % effective area [um^2]
n2        = 4.8e-20;                    % nonlinear index [m^2/W]
lambda    = link.lambda;                % wavelength [nm] @ dispersion
D         = -100;                       % dispersion [ps/nm/km] @ wavelength
S         = 0;                          % slope [ps/nm^2/km] @ wavelength
comp_LL   = -link.disp*link.LL/1e3/D*1e3;    % comp. fiber length [m];           

Ns_prop   = 100;                        % number of SSFM propagation step

comp_ch   = Channel(comp_LL,comp_alphadB,lambda,aeff,n2,D,S,Ns_prop);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                         DSP parameters                                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Ns_bprop      = sp.bprop;              % SSFM and ESSFM backpropagation steps
Ns_bprop_comp = 10;
dsp           = DSP(ch,Ns_bprop);
dsp_comp      = DSP(comp_ch,Ns_bprop_comp);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                      Signal parameters (wdm)                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
symbrate = signal.symbrate;           % symbol rate             [Gbaud]
Ps_dBm   = pdbm;                      % Power vector            [dBm]
Pavg     = 10.^(0.1*(Ps_dBm-30));     % Power vector            [W]
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
Gerbio01    = alphadB*LL*1e-3;
etasp01     = amp.etasp01;
amptype   = amp.type;
Gerbio02  = comp_alphadB*comp_LL*1e-3;
etasp02   = amp.etasp02;
amptype   = amp.type;
% ampli     = Ampliflat(Pavg,ch,Gerbio,etasp);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                         Matched filter                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Hf        = filt(pls,sig.FN);

set(sig,'POWER',Pavg);

Pu1dB = pdbm-4;
Pu1   = 10.^(0.1*(Pu1dB-30));
ampli01   = Ampliflat(Pu1,ch,Gerbio01,etasp01,amptype,Nspan);
ampli02   = Ampliflat(Pavg,ch,Gerbio02,etasp02,amptype,Nspan);

%% Subcarrier creation
for ii = 1:sig.NCH
    set(sub_sig{ii},'POWER',Pavg);
    Laser.GetLaserSource(Pavg,sub_sig{ii},lambda,0);
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
    
    sub_carriers{ii} = sub_sig{ii}.getproperties();
end

%% Sub carrier Multiplexing
if (not(sub_signal.nsc == 1))
    for ii = 1:sig.NCH
        MuxDemux.Mux(sub_sig{ii}.SUB_FIELDX,sub_sig{ii}.SUB_FIELDY,sub_sig{ii});
        sub_carriers{ii} = sub_sig{ii}.getproperties();
    end
end



%% Channel Multiplexing
% Laser.GetLaserSource(Pavg,sig,lambda,0.400835);
Laser.GetLaserSource(Pavg,sig,lambda,0);
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
%% Propagation
%     if gpu
%         gpu_propagation(ch,Nspan,ampli,sig);
%     else
%         propagation(ch,Nspan,ampli,sig);
%     end
    
   set(sig,'FIELDX', gpuArray(complex(get(sig,'FIELDX'))));
   set(sig,'FIELDY', gpuArray(complex(get(sig,'FIELDY'))));

    for i = 1:Nspan
        set(sig,'POWER',Pavg);
        sing_span_propagation(ch,sig,gpu);
        AddNoise(ampli01,sig);
        set(sig,'POWER',Pu1);
        sing_span_propagation(comp_ch,sig,gpu);
        AddNoise(ampli02,sig);        
    end
    
    set(sig,'POWER',Pavg);
    set(sig,'FIELDX', gather(get(sig,'FIELDX')));
    set(sig,'FIELDY', gather(get(sig,'FIELDY')));
            
%     g = gpuDevice(1);
%     reset(g);
    
%     for i = 1:Nspan
%         sing_span_propagation(ch,sig,'true')
%     end
%     AddNoise(ampli,sig);

%% Channel Demultiplexing
oHf       = myfilter(oftype,sig.FN,obw,0);    
[zfieldx,zfieldy] = MuxDemux.Demux(sig,oHf,cch);
set(sig,'FIELDX',zfieldx(:,cch));
set(sig,'FIELDY',zfieldy(:,cch));

%% Backpropagation

 %backpropagation(dsp,Pavg*10^(-Gerbio*0.1),sig,Nspan,'ssfm');
 
    set(sig,'FIELDX', gpuArray(complex(get(sig,'FIELDX'))));
    set(sig,'FIELDY', gpuArray(complex(get(sig,'FIELDY'))));
 
    for i = 1:Nspan
        backpropagation(dsp_comp,Pu1*10^(-Gerbio02*0.1),sig,1,'ssfm',0);
        backpropagation(dsp,Pavg*10^(-Gerbio01*0.1),sig,1,'ssfm',0);        
    end
    
    set(sig,'FIELDX', gather(get(sig,'FIELDX')));
    set(sig,'FIELDY', gather(get(sig,'FIELDY')));
    
    g      = gpuDevice(1);
    reset(g);

%     if (strcmp(amptype,'Raman'))
%         dsp.backpropagation(Pavg,sig,Nspan,'ssfm',gpu);
%     else
%         dsp.backpropagation(Pavg*10^(-Gerbio*0.1),sig,Nspan,'ssfm',gpu);
%     end

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
     Laser.GetLaserSource(Pavg,sig,lambda,0);
     if(pol == 2)
         [sig.SUB_FIELDX,sig.SUB_FIELDY] = MuxDemux.Demux(sig,oHf,0);
     else
         [sig.SUB_FIELDX] = MuxDemux.Demux(sig,oHf,0);
     end
     dsp.scdownsampling(sig);
 end
 
 
 
 signals      = sig.getproperties();
 Gm1          = (10^(Gerbio01*0.1)-1.d0)/Pu1 + (10^(Gerbio02*0.1)-1.d0)/Pavg ;
 N0           = etasp01*(Gm1)*HPLANCK*CLIGHT/(ch.lambda * 1e-9);
 SNRdB        = 10*log10(1/(symbrate*sub_signal.nsc)/10^9/N0/Nspan);
 
end

