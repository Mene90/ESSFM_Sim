function [ data ] = BER_ESSFM_XY(Nstep,NC,dBm,sym_length,n_prop_steps,etasp)
%% Example of FIELD propagation with single polarization
% This example calculate the ber of the received field after
% backpropagation. The programm is basically divided in slots. 
% In the first one the the Link, Tx and Rx parameters are seted


%% Initialization of channel, trx and rx
% Description of first code block

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                         Link parameters                                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
LL        = 1e5;                % length [m]
alphadB   = 0.2;                  % attenuation [dB/km]
aeff      = 80;                   % effective area [um^2]
n2        = 2.5e-20;              % nonlinear index [m^2/W]
lambda    = 1550;                 % wavelength [nm] @ dispersion
D         = 17;                   % dispersion [ps/nm/km] @ wavelength
S         = 0;                    % slope [ps/nm^2/km] @ wavelength

Ns_prop   = n_prop_steps;         % number of SSFM propagation step
Nspan     = 40;                   % total number of amplifiers

pmd       = false;                % pmd enable/disable

ch        = Channel(LL,alphadB,lambda,aeff,n2,D,S,Ns_prop,pmd);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                         DSP parameters                                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Ns_bprop = Nstep;                  % SSFM and ESSFM backpropagation steps


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                      Global Signal parameters                          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
symbrate  = 32;                  % symbol rate [Gbaud]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                         Pulse parameters                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Ps_dBm   = dBm;                    % Power vector            [dBm]
Pavg     = 10.^(0.1*(Ps_dBm-30));     % total transmitted power [W]
Plen     = length(Ps_dBm);

pls.shape   = 'RRC';
pls.bw      = 1.0;                       % duty cycle
pls.ord     = 0.2;                       % pulse roll-off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                         Amplifier parameters                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Gerbio    = alphadB*LL*1e-3;
% etasp     = 2;

%% Training phase
% The parameters for the ESSFM are calculating through a training signal
% with 2^10 symbols with a qpsk modulation

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                      Trainin Signal parameters                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Nsymb     = 2^10;                % number of symbols
Nt        = 2;                   % points x symbol
nfft      = Nsymb * Nt;
sig       = Signal(Nsymb,Nt,symbrate);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                         Matched filter                                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Hf       = gpuArray(transpose(filt(pls,sig.FN)));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                         ESSFM PARAMETERS                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

options = optimset('Algorithm','trust-region-reflective','Display','off',...
    'Jacobian','off','DerivativeCheck','off','TolFun',1e-13,'TolX',1e-13);

C0        = zeros(NC,1);
C0(1,1)   = 1;
Loss      = 10^(-Gerbio*0.1);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for nn=1:Plen
      
    dsp       = DSP(ch,Ns_bprop);
    sig       = Signal(Nsymb,Nt,symbrate);
    
    
    ampli     = Ampliflat(Pavg(nn),ch,Gerbio,etasp);
    
    [patx{nn}(:,1), patmatx]    = Pattern.debruijn(1,4,Nsymb);
    [paty{nn}(:,1), patmaty]    = Pattern.debruijn(2,4,Nsymb);
    
    E                  = Laser.GetLaserSource(Pavg(nn), nfft);
    
    set(sig,'POWER'     ,Pavg(nn));
    set(sig,'FIELDX_TX',Modulator.ApplyModulation(E, 2*patmatx-1, sig, pls));
    set(sig,'FIELDX'   ,Modulator.ApplyModulation(E, 2*patmatx-1, sig, pls));
    set(sig,'FIELDY_TX',Modulator.ApplyModulation(E, 2*patmaty-1, sig, pls));
    set(sig,'FIELDY'   ,Modulator.ApplyModulation(E, 2*patmaty-1, sig, pls));
    
    ux = gpuArray(complex(get(sig,'FIELDX')));
    uy = gpuArray(complex(get(sig,'FIELDY')));
    
    for i = 1:Nspan
        [ux, uy]      = ch.gpu_vectorial_ssfm(Pavg(nn),sig,ux,uy);
        [ux, uy]      = ampli.gpu_AddNoise(sig,ux,uy);
    end
    
     fmin      = @(C) gpu_vec_essfm_opt(sig,ux,uy,dsp,C,Nspan,Loss,Hf);

    [Coeff(nn,:),err]=lsqnonlin(fmin,C0,[],[],options);
end

%% Transmission Propagation and Reception
% Propagation and backpropagation of a signal of lenght 2^22 of random 
% symbols with qpsk modulation. After the propagation the BER is estimated
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           GPU optimization                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C = gpuArray(Coeff);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                      Trainin Signal parameters                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Nsymb     = sym_length;                % number of symbols
Nt        = 2;                         % points x symbol
nfft      = Nsymb * Nt;

sig       = Signal(Nsymb,Nt,symbrate);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                         Matched filter                                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Hf_BER        = transpose(filt(pls,sig.FN));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  
for nn=1:Plen
 
 
    dsp       = DSP(ch,Ns_bprop);
    sig       = Signal(Nsymb,Nt,symbrate);
        
    ampli                 = Ampliflat(Pavg(nn),ch,Gerbio,etasp);
    
    [patx_tx{nn}(:,1), patmatx_tx]    = Pattern.random(4,Nsymb);
    [paty_tx{nn}(:,1), patmaty_tx]    = Pattern.random(4,Nsymb);
    
    E                     = Laser.GetLaserSource(Pavg(nn), nfft);
        
    set(sig,'POWER'     ,Pavg(nn));
    set(sig,'FIELDX'    ,Modulator.ApplyModulation(E, 2*patmatx_tx-1, sig, pls));
    set(sig,'FIELDX_TX' ,Modulator.ApplyModulation(E, 2*patmatx_tx-1, sig, pls));
    set(sig,'FIELDY'    ,Modulator.ApplyModulation(E, 2*patmaty_tx-1, sig, pls));
    set(sig,'FIELDY_TX' ,Modulator.ApplyModulation(E, 2*patmaty_tx-1, sig, pls));
    
    
    ux = gpuArray(complex(get(sig,'FIELDX')));
    uy = gpuArray(complex(get(sig,'FIELDY')));
    for i = 1:Nspan
        [ux, uy]      = ampli.gpu_AddNoise(sig,ux,uy);
        [ux, uy]      = ch.gpu_vectorial_ssfm(Pavg(nn),sig,ux,uy);
    end

    ux_enh = ux;
    uy_enh = uy;
    
    sig_st_rx = copy(sig);
    sig_enh_rx = copy(sig);
    if(Nstep>=1)
        for i = 1:Nspan
            [ux, uy]           = dsp.DBP_gpu_vec_ssfm (Pavg(nn)*Loss,sig_st_rx,ux,uy);
            [ux_enh, uy_enh]   = dsp.DBP_gpu_vec_essfm(Pavg(nn)*Loss,sig_enh_rx,C(nn,:)',ux_enh,uy_enh);
        end
    else
        for i = 1:round(Nspan*Nstep)
            [ux, uy]           = dsp.DBP_gpu_vec_ssfm (Pavg(nn)*Loss,sig_st_rx,ux,uy);
            [ux_enh, uy_enh]   = dsp.DBP_gpu_vec_essfm(Pavg(nn)*Loss,sig_enh_rx,C(nn,:)',ux_enh,uy_enh);
        end
    end
    
    set(sig_st_rx,'FIELDX',gather(ux));
    set(sig_st_rx,'FIELDY',gather(uy));
    set(sig_enh_rx,'FIELDX',gather(ux_enh));
    set(sig_enh_rx,'FIELDY',gather(uy_enh));
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %                           X-Pol                                    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    FIELDX_TX       = get(sig       ,'FIELDX_TX');  
    FIELDX_ST_RX    = get(sig_st_rx ,'FIELDX'   );
    FIELDX_ENH_RX   = get(sig_enh_rx,'FIELDX'   );
    
    FIELDX_ST_RX    = ifft(fft(FIELDX_ST_RX).*Hf_BER);
    FIELDX_ENH_RX   = ifft(fft(FIELDX_ENH_RX).*Hf_BER);
    
    rotx_st         = angle(mean(FIELDX_ST_RX .*conj(FIELDX_TX)));
    rotx_enh        = angle(mean(FIELDX_ENH_RX.*conj(FIELDX_TX)));
    
    FIELDX_ST_RX    = FIELDX_ST_RX *exp(-1i*rotx_st);
    FIELDX_ENH_RX   = FIELDX_ENH_RX*exp(-1i*rotx_enh);
    
    patmatx_st_rx    = samp2pat(angle(FIELDX_ST_RX(1:Nt:end)));
    patmatx_enh_rx   = samp2pat(angle(FIELDX_ENH_RX(1:Nt:end)));
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %                           Y-Pol                                    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    FIELDY_TX       = get(sig       ,'FIELDY_TX');
    FIELDY_ST_RX    = get(sig_st_rx ,'FIELDY'   );
    FIELDY_ENH_RX   = get(sig_enh_rx,'FIELDY'   );
    
    FIELDY_ST_RX    = ifft(fft(FIELDY_ST_RX).*Hf_BER);
    FIELDY_ENH_RX   = ifft(fft(FIELDY_ENH_RX).*Hf_BER);    
    
    roty_st         = angle(mean(FIELDY_ST_RX .*conj(FIELDY_TX)));
    roty_enh        = angle(mean(FIELDY_ENH_RX.*conj(FIELDY_TX))); 
    
    FIELDY_ST_RX    = FIELDY_ST_RX *exp(-1i*roty_st);
    FIELDY_ENH_RX   = FIELDY_ENH_RX*exp(-1i*roty_enh);    
   
    patmaty_st_rx    = samp2pat(angle(FIELDY_ST_RX(1:Nt:end)));
    patmaty_enh_rx   = samp2pat(angle(FIELDY_ENH_RX(1:Nt:end)));
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    avgberx = [ber(patmatx_st_rx, patmatx_tx) ber(patmatx_enh_rx,patmatx_tx)];
    avgbery = [ber(patmaty_st_rx, patmaty_tx) ber(patmaty_enh_rx,patmaty_tx)];
    
    avgber = 0.5*(avgberx+avgbery);
    
    data(nn,:) = avgber;
    
%     display(['Power[dBm] = '   , num2str(dBm(nn))  ,char(9)...
%              'BER con SSFM = ' , num2str(avgber(1)),char(9)...
%              'BER con ESSFM = ', num2str(avgber(2))            ]);

end


end  