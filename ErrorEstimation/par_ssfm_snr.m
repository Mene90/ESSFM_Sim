function [ avg_snr ] = par_ssfm_snr( Ps_dBm,ch,dsp,SNR_sig,ampli_par,system)
sig   = copy(SNR_sig);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                         Some System Parameters                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Pavg      = 10.^(0.1*(Ps_dBm-30));            % total transmitted power [W]
Nspan     = system.Nspan;
Loss      = system.Loss;
Hf_SNR    = system.mfil_snr;
set(sig  ,'POWER'     ,Pavg);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                         Amplifier parameters                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ampli     = Ampliflat(Pavg,ch,ampli_par.G,ampli_par.e);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                      PROPAGATION                                        % 
%                    + BACKPROPAGATION SSFM                               % 
%                    + SNR ESTIMATION                                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i = 1:Nspan
    sig      = ch.vectorial_ssfm(Pavg,sig);
    sig      = ampli.AddNoise(sig);
end

sig_st_rx = copy(sig);
for i = 1:Nspan
    sig_st_rx   = dsp.DBP_vec_ssfm(Pavg*Loss,sig_st_rx);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           X-Pol                                    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

FIELDX_TX       = get(sig       ,'FIELDX_TX');
FIELDX_ST_RX    = get(sig_st_rx ,'FIELDX'   );

FIELDX_ST_RX    = ifft(fft(FIELDX_ST_RX).*Hf_SNR);

rotx_enh        = angle(mean(FIELDX_ST_RX.*conj(FIELDX_TX)));

FIELDX_ST_RX    = FIELDX_ST_RX*exp(-1i*rotx_enh);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           Y-Pol                                    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

FIELDY_TX       = get(sig       ,'FIELDY_TX');
FIELDY_ST_RX    = get(sig_st_rx ,'FIELDY'   );

FIELDY_ST_RX    = ifft(fft(FIELDY_ST_RX).*Hf_SNR);

roty_enh        = angle(mean(FIELDY_ST_RX.*conj(FIELDY_TX)));

FIELDY_ST_RX    = FIELDY_ST_RX*exp(-1i*roty_enh);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

snr_x = SNR(FIELDX_ST_RX(1:sig.NT:end),FIELDX_TX(1:sig.NT:end));
snr_y = SNR(FIELDY_ST_RX(1:sig.NT:end),FIELDY_TX(1:sig.NT:end));

avg_snr = 10*log10(0.5*(snr_x+snr_y));

end
