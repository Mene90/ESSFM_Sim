function [ data ] = BER_ESSFM_X(Nstep,NC,dBm,sym_length,n_prop_steps,etasp)
%% Example of FIELD propagation with single polarization
% This example calculate the ber of the received field after
% backpropagation. The programm is basically divided in slots. 
% In the first one the the Link, Tx and Rx parameters are seted

%% Initialization of channel, tx and rx
% In this first block the link, amplifier, DSP, signal and pulse
% parameters are initialized  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                         Link parameters                                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
LL        = 1.2e5;                % length [m]
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

Ns_bprop = Nstep;                   % SSFM and ESSFM backpropagation steps


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                      Global Signal parameters                          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
symbrate  = 32;                       % symbol rate [Gbaud]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                         Pulse parameters                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Ps_dBm   = dBm;                       % Power vector            [dBm]
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

Hf       = transpose(filt(pls,sig.FN));

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
    
    dsp     = DSP(ch,Ns_bprop);    
    sig       = Signal(Nsymb,Nt,symbrate);
    fmin      = @(C) essfm_opt(sig,dsp,C,Nspan,Loss,Hf);
    
    ampli                 = Ampliflat(Pavg(nn),ch,Gerbio,etasp);
    
    [pat{nn}(:,1), patmat]    = Pattern.debruijn(1,4,Nsymb);
    
    E                     = Laser.GetLaserSource(Pavg(nn), nfft);
    
    set(sig,'POWER'     ,Pavg(nn));
    set(sig,'FIELDX_TX',Modulator.ApplyModulation(E, 2*patmat-1, sig, pls));
    set(sig,'FIELDX'   ,Modulator.ApplyModulation(E, 2*patmat-1, sig, pls));
    
    for i = 1:Nspan
        sig      = ampli.AddNoise(sig);
        sig      = ch.scalar_ssfm(Pavg(nn),sig);
        
    end
    
    [C(nn,:),err] = lsqnonlin(fmin,C0,[],[],options);
end

%% Transmission Propagation and Reception
% Propagation and backpropagation of a signal of lenght 2^22 of random 
% symbols with qpsk modulation. After the propagation the BER is estimated

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           Signal parameters                            %
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
    
    dsp     = DSP(ch,Ns_bprop);
    sig       = Signal(Nsymb,Nt,symbrate);
    
    ampli                     = Ampliflat(Pavg(nn),ch,Gerbio,etasp);
    
    [patx{nn}(:,1), patmat_tx]    = Pattern.random(4,Nsymb);
    
    E                         = Laser.GetLaserSource(Pavg(nn), nfft);
    
    set(sig,'POWER'     ,Pavg(nn));
    set(sig,'FIELDX'    ,Modulator.ApplyModulation(E, 2*patmat_tx-1, sig, pls));
    set(sig,'FIELDX_TX' ,Modulator.ApplyModulation(E, 2*patmat_tx-1, sig, pls));
    tic
    for i = 1:Nspan
         sig      = ampli.AddNoise(sig);
         sig      = ch.scalar_ssfm(Pavg(nn),sig);     
    end
    toc
    set(sig,'FIELDX'    ,sig.FIELDX_TX);
    
    tic
    ux = gpuArray(complex(get(sig,'FIELDX')));
    for i = 1:Nspan
         ux      = ampli.gpu_AddNoise(sig,ux);
         ux      = ch.gpu_scalar_ssfm(Pavg(nn),sig,ux);     
    end
    set(sig,'FIELDX',gather(ux));
    toc
    
    
     sig_st_rx = copy(sig);
    sig_enh_rx = copy(sig);
    
    if (Nstep>=1)
        for i = 1:Nspan
            sig_st_rx    = dsp.DBP_scalar_ssfm(Pavg(nn)*Loss,sig_st_rx);
            sig_enh_rx  = dsp.DBP_essfm(Pavg(nn)*Loss,sig_enh_rx,C(nn,:)');
        end
    else
        for i = 1:round(Nspan*Nstep)
            sig_st_rx    = dsp.DBP_scalar_ssfm(Pavg(nn)*Loss,sig_st_rx);
            sig_enh_rx  = dsp.DBP_essfm(Pavg(nn)*Loss,sig_enh_rx,C(nn,:)');
        end
    end
    
    
    
    FIELDX_TX       = get(sig,'FIELDX_TX');
    FIELDX_ST_RX    = get(sig_st_rx,'FIELDX');
    FIELDX_ENH_RX   = get(sig_enh_rx,'FIELDX');
    
    
    FIELDX_ST_RX    = ifft(fft(FIELDX_ST_RX).*Hf_BER);
    FIELDX_ENH_RX   = ifft(fft(FIELDX_ENH_RX).*Hf_BER);
    
    rot_st          = angle(mean(FIELDX_ST_RX .*conj(FIELDX_TX)));
    rot_enh         = angle(mean(FIELDX_ENH_RX.*conj(FIELDX_TX)));
    FIELDX_ST_RX    = FIELDX_ST_RX *exp(-1i*rot_st);
    FIELDX_ENH_RX   = FIELDX_ENH_RX*exp(-1i*rot_enh);
    
    patmat_st_rx    = samp2pat(angle(FIELDX_ST_RX(1:Nt:end)));
    patmat_enh_rx   = samp2pat(angle(FIELDX_ENH_RX(1:Nt:end)));
    
    avgber = [ber(patmat_st_rx,patmat_tx) ber(patmat_enh_rx,patmat_tx)];
    
    data(nn,:) = avgber;
    display(['Power[dBm] = '   , num2str(dBm(nn))  ,char(9)...
             'BER con SSFM = ' , num2str(avgber(1)),char(9)...
             'BER con ESSFM = ', num2str(avgber(2))            ]);
end


end  