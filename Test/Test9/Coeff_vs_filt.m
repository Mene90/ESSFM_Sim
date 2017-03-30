function [ data ] = Coeff_vs_filt( Nstep,NC,dBm,sym_length,n_prop_steps,etasp )
%% Example of FIELD propagation with single polarization
% This example calculate the ber of the received field after
% backpropagation. The programm is basically divided in slots. 
% In the first one the the Link, Tx and Rx parameters are seted


%% Initialization of channel, trx and rx
% Description of first code block

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                         Link parameters                                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
LL        = 1e5;                  % length [m]
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
symbrate  = 50;                  % symbol rate [Gbaud]

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
Nsymb     = 2^16;                % number of symbols
Nt        = 2;                   % points x symbol
nfft      = Nsymb * Nt;
sig       = Signal(Nsymb,Nt,symbrate);
trainlength = Nsymb;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                         Matched filter                                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Hf       = gpuArray(transpose(filt(pls,sig.FN)));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                         ESSFM PARAMETERS                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

C0        = zeros(NC,1);
C0(1,1)   = 1;
Loss      = 10^(-Gerbio*0.1);

options = optimset('Algorithm','trust-region-reflective','Display','iter',...
    'Jacobian','off','DerivativeCheck','off','TolFun',1e-13,'TolX',1e-13,'MaxFunEvals',10000*length(C0));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
bw = 0.02;
seedx = 1;
seedy = sig.NSYMB/2^3;      
    dsp       = DSP(ch,Ns_bprop);
    sig       = Signal(Nsymb,Nt,symbrate);
    
    
    ampli     = Ampliflat(Pavg(1),ch,Gerbio,etasp);
    
    [patx{1}(:,1), patmatx]    = Pattern.debruijn(seedx,4,Nsymb);
    [paty{1}(:,1), patmaty]    = Pattern.debruijn(seedy,4,Nsymb);
    
    E                  = Laser.GetLaserSource(Pavg(1), nfft);
    
    set(sig,'POWER'     ,Pavg(1));
    set(sig,'FIELDX_TX',1./sqrt(2.)*((2*patmatx(:,1)-1)+1i*(2.*patmatx(:,2)-1)));
    set(sig,'FIELDY_TX',1./sqrt(2.)*((2*patmaty(:,1)-1)+1i*(2.*patmaty(:,2)-1)));
    set(sig,'FIELDX'   ,Modulator.ApplyModulation(E, 2*patmatx-1, sig, pls));
    set(sig,'FIELDY'   ,Modulator.ApplyModulation(E, 2*patmaty-1, sig, pls));
    
    ux = gpuArray(complex(get(sig,'FIELDX')));
    uy = gpuArray(complex(get(sig,'FIELDY')));
    
    for i = 1:Nspan
        [ux, uy]      = ch.gpu_vectorial_ssfm(Pavg(1),sig,ux,uy);
        [ux, uy]      = ampli.gpu_AddNoise(sig,ux,uy);
    end
    
     fmin      = @(C) gpu_vec_essfm_opt(sig,ux,uy,dsp,C,Nspan,Loss,Hf);

    [Coeff,err]=lsqnonlin(fmin,C0,[],[],options);

Coeff.'
%% Transmission Propagation and Reception
% Propagation and backpropagation of a signal of lenght 2^22 of random 
% symbols with qpsk modulation. After the propagation the BER is estimated
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           GPU optimization                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% C = gpuArray(Coeff);
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

% Hf_BER        = transpose(filt(pls,sig.FN));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pls.shape = 'gauss';
pls.ord   = 2;
pls.bw    = bw;
H = ifft(transpose(filt(pls,sig.FN)));
Hf = [H(length(H)-(length(Coeff)-2):end);H(1:length(Coeff))];
CC = [flipud(Coeff);Coeff(2:end)];
fig=figure(1);

stem([(-(length(Coeff)-1):length(Coeff)-1)',(-(length(Coeff)-1):length(Coeff)-1)'],[Hf,CC]);
grid on;
lg = {strcat('ESSFM Coeff =',{' '},int2str(length(Coeff)-1))};
legend('Gauss Filt',lg{1,1}{1,1});

savefig(fig,strcat('plot/Ns',...
    int2str(dsp.nstep*Nspan),...
    '_Nc',int2str(length(Coeff)-1),'_dBm',int2str(Ps_dBm(1)),'_lengthtrain',int2str(log2(trainlength)),'.fig'));
hold('off');
close(fig);

end

