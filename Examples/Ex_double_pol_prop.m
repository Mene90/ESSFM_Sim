function [ data ] = Ex_sing_pol_prop(Nstep,NC)
%% Example of FIELD propagation with single polarization
% This example calculate the ber of the received field after
% backpropagation. The programm is basically divided in slots. 
% In the first one the the Link, Tx and Rx parameters are seted

addpath('C:\Users\mene9\Documents\MATLAB\My Simulator\Link')
addpath('C:\Users\mene9\Documents\MATLAB\My Simulator\Rx')
addpath('C:\Users\mene9\Documents\MATLAB\My Simulator\Tx')
addpath('C:\Users\mene9\Documents\MATLAB\My Simulator\ErrorEstimation')
addpath('C:\Users\mene9\Documents\MATLAB\My Simulator\')

%% Initialization of channel, trx and rx
% Description of first code block

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                         Link parameters                                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
LL        = 1.2e5;                % length [m]
alphadB   = 0.2;                  % attenuation [dB/km]
aeff      = 87;                   % effective area [um^2]
n2        = 2.5e-20;              % nonlinear index [m^2/W]
lambda    = 1550;                 % wavelength [nm] @ dispersion
D         = 17;                   % dispersion [ps/nm/km] @ wavelength
S         = 0;                    % slope [ps/nm^2/km] @ wavelength

Ns_prop   = 10;                   % number of SSFM propagation step
Nspan     = 40;                   % total number of amplifiers

pmd       = false;                % pmd enable/disable

ch        = Channel(LL,alphadB,lambda,aeff,n2,D,S,Ns_prop,pmd);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                         DSP parameters                                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Ns_bprop = Nstep;                  % SSFM and ESSFM backpropagation steps

dsp     = DSP(ch,Ns_bprop);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                      Global Signal parameters                          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
symbrate  = 32;                  % symbol rate [Gbaud]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                         Pulse parameters                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Ps_dBm   = (-2:4);                    % Power vector            [dBm]
Pavg     = 10.^(0.1*(Ps_dBm-30));     % total transmitted power [W]
Plen     = length(Ps_dBm);

duty     = 1.0;                       % duty cycle
roll     = 0.2;                       % pulse roll-off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                         Amplifier parameters                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Gerbio    = alphadB*LL*1e-3;
etasp     = 2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                         Matched filter                                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Training phase
% Description of second code block

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                      Trainin Signal parameters                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Nsymb     = 2^10;                % number of symbols
Nt        = 2;                   % points x symbol
nfft      = Nsymb * Nt;

sig       = Signal(Nsymb,Nt,symbrate);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                         ESSFM PARAMETERS                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

options = optimset('Algorithm','trust-region-reflective','Display','off',...
    'Jacobian','off','DerivativeCheck','off','TolFun',1e-13,'TolX',1e-13);

C0        = zeros(NC,1);
C0(1,1)   = 1;
Loss      = 10^(-Gerbio*0.1);
fmin      = @(C) vec_essfm_opt(sig,dsp,C,Nspan,Loss);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for nn=1:Plen
    
    ampli     = Ampliflat(Pavg(nn),ch,Gerbio,etasp);
    
    [patx(:,1), patmatx]    = Pattern.debruijn(1,4,Nsymb);
    [paty(:,1), patmaty]    = Pattern.debruijn(2,4,Nsymb);
    
    E                  = Laser.GetLaserSource(Pavg(nn), nfft);
    v                  = ElectricSource(sig,'cosroll',duty,roll);
    
    elecx_i            = v.pat2electricsource(patmatx(:,1),'qpsk');
    elecx_q            = v.pat2electricsource(patmatx(:,2),'qpsk');
    elecy_i            = v.pat2electricsource(patmaty(:,1),'qpsk');
    elecy_q            = v.pat2electricsource(patmaty(:,2),'qpsk');
    
    set(sig,'POWER'     ,Pavg(nn));
    set(sig,'FIELDX_TX',QI_modulator.ApplyModulation(E, elecx_i, elecx_q));
    set(sig,'FIELDX'   ,QI_modulator.ApplyModulation(E, elecx_i, elecx_q));
    set(sig,'FIELDY_TX',QI_modulator.ApplyModulation(E, elecy_i, elecy_q));
    set(sig,'FIELDY'   ,QI_modulator.ApplyModulation(E, elecy_i, elecy_q));
    
    for i = 1:Nspan
        sig      = ch.vectorial_ssfm(Pavg(nn),sig);
        sig      = ampli.AddNoise(sig);
    end
    
    [C(nn,:),err]=lsqnonlin(fmin,C0,[],[],options);
end
%% Transmission Propagation and Reception
% Description of second code block

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                      Trainin Signal parameters                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Nsymb     = 2^16;                % number of symbols
Nt        = 2;                   % points x symbol
nfft      = Nsymb * Nt;

sig       = Signal(Nsymb,Nt,symbrate);

for nn=1:Plen
    
    ampli                     = Ampliflat(Pavg(nn),ch,Gerbio,etasp);
    
    [patx_tx(:,1), patmatx_tx]    = Pattern.random(4,Nsymb);
    [paty_tx(:,1), patmaty_tx]    = Pattern.random(4,Nsymb);
    
    E                     = Laser.GetLaserSource(Pavg(nn), nfft);
    v                     = ElectricSource(sig,'cosroll',duty,roll);
    
    elecx_i               = v.pat2electricsource(patmatx_tx(:,1),'qpsk');
    elecx_q               = v.pat2electricsource(patmatx_tx(:,2),'qpsk');
    elecy_i               = v.pat2electricsource(patmaty_tx(:,1),'qpsk');
    elecy_q               = v.pat2electricsource(patmaty_tx(:,2),'qpsk');
    
    set(sig,'POWER'     ,Pavg(nn));
    set(sig,'FIELDX'    ,QI_modulator.ApplyModulation(E, elecx_i, elecx_q));
    set(sig,'FIELDX_TX' ,QI_modulator.ApplyModulation(E, elecx_i, elecx_q));    
    set(sig,'FIELDY_TX' ,QI_modulator.ApplyModulation(E, elecy_i, elecy_q));
    set(sig,'FIELDY'    ,QI_modulator.ApplyModulation(E, elecy_i, elecy_q));
    
    
    for i = 1:Nspan
        sig      = ch.vectorial_ssfm(Pavg(nn),sig);
        sig      = ampli.AddNoise(sig);
    end
    
    sig_st_rx = copy(sig);
    for i = 1:Nspan
        sig_st_rx    = dsp.DBP_vec_ssfm(Pavg(nn)*Loss,sig_st_rx);
    end
    
    sig_enh_rx = copy(sig);
    for i = 1:Nspan
        sig_enh_rx  = dsp.DBP_vec_essfm(Pavg(nn)*Loss,sig_enh_rx,C(nn));
    end
    
    FIELDX_TX       = get(sig       ,'FIELDX_TX');
    FIELDY_TX       = get(sig       ,'FIELDY_TX');
    FIELDX_ST_RX    = get(sig_st_rx ,'FIELDX'   );
    FIELDX_ENH_RX   = get(sig_enh_rx,'FIELDX'   );
    FIELDY_ST_RX    = get(sig_st_rx ,'FIELDY'   );
    FIELDY_ENH_RX   = get(sig_enh_rx,'FIELDY'   );
    
    rotx_st         = angle(mean(FIELDX_ST_RX .*conj(FIELDX_TX)));
    roty_st         = angle(mean(FIELDY_ST_RX .*conj(FIELDY_TX)));
    rotx_enh        = angle(mean(FIELDX_ENH_RX.*conj(FIELDX_TX)));
    roty_enh        = angle(mean(FIELDY_ENH_RX.*conj(FIELDY_TX)));
    
    FIELDX_ST_RX    = FIELDX_ST_RX *exp(-1i*rotx_st);
    FIELDY_ST_RX    = FIELDY_ST_RX *exp(-1i*roty_st);
    FIELDX_ENH_RX   = FIELDX_ENH_RX*exp(-1i*rotx_enh);
    FIELDY_ENH_RX   = FIELDY_ENH_RX*exp(-1i*roty_enh);
    
    patmatx_st_rx    = samp2pat(angle(FIELDX_ST_RX(1:Nt:end)));
    patmaty_st_rx    = samp2pat(angle(FIELDY_ST_RX(1:Nt:end)));
    
    patmatx_enh_rx   = samp2pat(angle(FIELDX_ENH_RX(1:Nt:end)));
    patmaty_enh_rx   = samp2pat(angle(FIELDY_ENH_RX(1:Nt:end)));
    
    avgberx = [ber(patmatx_st_rx, patmatx_tx) ber(patmatx_enh_rx,patmatx_tx)];
    avgbery = [ber(patmaty_st_rx, patmaty_tx) ber(patmaty_enh_rx,patmaty_tx)];
    
    avgber = 0.5*(avgberx+avgbery);
    
    display(['BER con SSFM = ', num2str(avgber(1)) ,char(9) 'BER con ESSFM = ',num2str(avgber(2))]);
    
end

clear;

end  