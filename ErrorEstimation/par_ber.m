function [ avgber ] = par_ber(Ps_dBm,ch,dsp,ber_sig,tr_sig,pat_tx,ampli_par,system,essfm)
sig   = copy(ber_sig);
t_sig = copy(tr_sig);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                         Some System Parameters                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Pavg      = 10.^(0.1*(Ps_dBm-30));            % total transmitted power [W]
Nspan     = system.Nspan;
Loss      = system.Loss;
Hf_BER    = system.mfil_ber;
Hf        = system.mfil;
set(t_sig,'POWER'     ,Pavg);
set(sig  ,'POWER'     ,Pavg);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                         Amplifier parameters                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ampli     = Ampliflat(Pavg,ch,ampli_par.G,ampli_par.e);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                            BER Parameters                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
patmatx_tx = pat_tx.x;
patmaty_tx = pat_tx.y;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                         ESSFM PARAMETERS                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fmin      = @(C) vec_essfm_opt(t_sig,dsp,C,Nspan,Loss,Hf);
options   = essfm.opt;
C0        = zeros(essfm.NC,1);
C0(1,1)   = 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                               TRAINING                                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    for i = 1:Nspan
        t_sig      = ch.vectorial_ssfm(Pavg,t_sig);
        t_sig      = ampli.AddNoise(t_sig);
    end
    
    [C(1,:),err]=lsqnonlin(fmin,C0,[],[],options);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                             BER ESTIMATION                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i = 1:Nspan
    sig      = ch.vectorial_ssfm(Pavg,sig);
    sig      = ampli.AddNoise(sig);
end

sig_enh_rx = copy(sig);
for i = 1:Nspan
    sig_enh_rx   = dsp.DBP_vec_essfm(Pavg*Loss,sig_enh_rx,C(1,:)');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           X-Pol                                    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

FIELDX_TX       = get(sig       ,'FIELDX_TX');
FIELDX_ENH_RX   = get(sig_enh_rx,'FIELDX'   );

FIELDX_ENH_RX   = ifft(fft(FIELDX_ENH_RX).*Hf_BER);

rotx_enh        = angle(mean(FIELDX_ENH_RX.*conj(FIELDX_TX)));

FIELDX_ENH_RX   = FIELDX_ENH_RX*exp(-1i*rotx_enh);

patmatx_enh_rx   = samp2pat(angle(FIELDX_ENH_RX(1:sig.NT:end)));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           Y-Pol                                    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

FIELDY_TX       = get(sig       ,'FIELDY_TX');
FIELDY_ENH_RX   = get(sig_enh_rx,'FIELDY'   );

FIELDY_ENH_RX   = ifft(fft(FIELDY_ENH_RX).*Hf_BER);

roty_enh        = angle(mean(FIELDY_ENH_RX.*conj(FIELDY_TX)));

FIELDY_ENH_RX   = FIELDY_ENH_RX*exp(-1i*roty_enh);

patmaty_enh_rx   = samp2pat(angle(FIELDY_ENH_RX(1:sig.NT:end)));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


avgberx = ber(patmatx_enh_rx,patmatx_tx);
avgbery = ber(patmaty_enh_rx,patmaty_tx);

avgber = log(0.5*(avgberx+avgbery)+1e-9);
end

