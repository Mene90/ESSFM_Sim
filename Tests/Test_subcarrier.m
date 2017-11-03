function [signals,sub_carriers,SNRdB,ch] = Test_subcarrier(link,sp,signal,sub_signal,amp,pdbm,wdm,pls,pol,gpu)
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
Ns_bprop = sp.bprop;              % SSFM and ESSFM backpropagation steps
dsp      = DSP(ch,Ns_bprop);
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
    Laser.GetLaserSource(Pavg,sub_sig{ii},lambda,0);
%     if sub_signal.nsc == 4
%         Laser.GetLaserSource(Pavg,sub_sig{ii},lambda,0.10017);
%     elseif sub_signal.nsc == 2
%         Laser.GetLaserSource(Pavg,sub_sig{ii},lambda,0.4005);
%     else
%         Laser.GetLaserSource(Pavg,sub_sig{ii},lambda,0.0665);
%     end
    
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
%     if (sub_signal.nsc == 4 || sub_signal.nsc == 6)
%         for ii = 1:sig.NCH
%             MuxDemux.Mux(sub_sig{ii}.SUB_FIELDX,sub_sig{ii}.SUB_FIELDY,sub_sig{ii},2);
%             
%         end
%     else
%         for ii = 1:sig.NCH
%             MuxDemux.Mux(sub_sig{ii}.SUB_FIELDX,sub_sig{ii}.SUB_FIELDY,sub_sig{ii},0);
%         end
%     end
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

MuxDemux.Mux(Eoptx,Eopty,sig);

%% Propagation
    if gpu
        gpu_propagation(ch,Nspan,ampli,sig);
    else
        propagation(ch,Nspan,ampli,sig);
    end
    
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

    if (strcmp(amptype,'Raman'))
        dsp.backpropagation(Pavg,sig,Nspan,'ssfm',gpu);
    else
        dsp.backpropagation(Pavg*10^(-Gerbio*0.1),sig,Nspan,'ssfm',gpu);
    end

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
     [sig.SUB_FIELDX,sig.SUB_FIELDY] = MuxDemux.Demux(sig,oHf,0);
%      if sub_signal.nsc == 4
% %          Laser.GetLaserSource(Pavg,sig,lambda,0.10017);
%          Laser.GetLaserSource(Pavg,sig,lambda,0);
%          [sig.SUB_FIELDX,sig.SUB_FIELDY] = MuxDemux.Demux(sig,oHf,0,2);
%      else
% %          Laser.GetLaserSource(Pavg,sig,lambda,0.4005);
%          Laser.GetLaserSource(Pavg,sig,lambda,0);
%          [sig.SUB_FIELDX,sig.SUB_FIELDY] = MuxDemux.Demux(sig,oHf,0,0);
%      end
     
     dsp.scdownsampling(sig);
 end
 
 
 
 signals      = sig.getproperties();
 SNRdB        = 10*log10(1/(symbrate*sub_signal.nsc)/10^9/ampli.N0/Nspan);
 
end

