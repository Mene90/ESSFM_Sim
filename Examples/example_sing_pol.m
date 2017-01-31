%% Example Title
% Summary of example objective
addpath('C:\Users\mene9\Documents\MATLAB\My Simulator\RX')
addpath('.\TX')
%% Section 1 Title
% Description of first code block
    %%%%%%%%%%%%%%%% Link parameters
    LL  = 1.2e5;         % length [m]
    alphadB = 0.2;     % attenuation [dB/km]
    aeff    = 80;      % effective area [um^2]
    n2      = 2.6e-20; % nonlinear index [m^2/W]
    lambda  = 1550;    % wavelength [nm] @ dispersion 
    D       = 17;      % dispersion [ps/nm/km] @ wavelength
    S       = 0;       % slope [ps/nm^2/km] @ wavelength
    dz      = LL/10;
    pmd     = false;
    
    symbrate = 32;     % symbol rate [Gbaud]
    Nsymb = 2^6;  % number of symbols
    Nt    = 2;   % points x symbol
    nfft    = Nsymb * Nt;
    
    %%%%%%%%%%%%%%%%%  Pulse parameters
    Pavg     = 0.0006;      % total transmitted power [W]
    duty     = 1.0;         % duty cycle
    roll     = 0.2;         % pulse roll-off


    [patx(:,1), patmatx]        = Pattern.debruijn(1,4,Nsymb);
    
    E = Laser.GetLaserSource(Pavg, nfft);
    sig = Signal(Nsymb,Nt,symbrate);
    v = ElectricSource(sig,'rootcosroll',duty,roll);
    
    elecx_i = v.pat2electricsource(patmatx(:,1),'qpsk');
    elecx_q = v.pat2electricsource(patmatx(:,2),'qpsk');
    
    Eoptx   = QI_modulator.ApplyModulation(E, elecx_i, elecx_q);
      
    sig.FIELDX = Eoptx ;
    
    x = 1:nfft;
%     plot(x,Eoptx,'-or',x,sig.FIELDX,'-g');
        
    ch      = Channel(LL,alphadB,lambda,aeff,n2,Pavg,D,S,dz,sig,pmd);
     
    [sig.FIELDX] = ch.ssfm(Eoptx);

     plot(x,Eoptx,'r',x,sig.FIELDX,'-g');
    
%     G=exp(ch.alphalin*LL);
%     
%     sig = ampli(Pavg,ch,G,2);
%     
%     plot(x,Eoptx,'-or',x,sig.FIELDX,'-g');
%     Nfft = length(sig.FN);
%     LO_Detuning = 0;
%     LO_PhaseNoise = zeros( Nfft, 1);
%     LO_Ecw = 1;
%     
%     rx      = receiver_cohmix(LO_Ecw,LO_PhaseNoise,LO_Detuning);
%     
%     ofilt.type ='gauss';
%     ofilt.obw  =1.9;
%     ofilt.oord =0;
%     
%     efilt.type = 'bessel5';
%     efilt.ebw = 0.65;
%     efilt.eord =0;
%     
%     Irx     = rx.receive(ch,ofilt,efilt);

      dsp     = DSP(ch);

    
%   [sig.FIELDX,sig.FIELDY] = dsp.A2D(Irx,Nt);
    
%     C0=zeros(4,1);
%     C0(1,1)=1;%-gm*Leff*Ps(nn);
%     options = optimset('Algorithm','trust-region-reflective','Display','off','Jacobian','off','DerivativeCheck','off','TolFun',1e-13,'TolX',1e-13);
%     [C,err]=lsqnonlin(@(C)essfm_opt(sig,patmatx,patmaty,dsp,C),C0,[],[],options);
%     
%     [FIELDX,FIELDY] = dsp.DBP_essfm(sig.FIELDX,sig.FIELDY,C);
% 
%     subplot(2,1,1);
%     plot(x,Eoptx,'r',x,FIELDX,'-g');
%     
%     subplot(2,1,2);
%     plot(x,Eopty,'r',x,FIELDY,'-g');

    [FIELDX] = dsp.DBP_ssfm(sig.FIELDX);
     
    %subplot(2,1,1);
    plot(x,Eoptx,'-or',x,FIELDX,'-g');
    
    clear
    clc
%     
%     Phases      = angle( [sig.FIELDX sig.FIELDY] );
%     subplot(2,1,1)
%     plot(x,angle(Eoptx),'r');
%     subplot(2,1,2);
%     plot(y,angle(sig.FIELDX),'g');
%     create_field('unique',Eoptx,Eopty,struct('power','average'));
    
%    patmat_hat = samp2pat(Phases);


%% Section 2 Title
% Description of second code block

