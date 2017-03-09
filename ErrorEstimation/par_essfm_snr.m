function [ avg_snr ] = par_essfm_snr( Ps_dBm,ch,dsp,SNR_sig,tr_sig,ampli_par,system,essfm )
sig   = copy(SNR_sig);
t_sig = copy(tr_sig);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                         Some System Parameters                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Pavg      = 10.^(0.1*(Ps_dBm-30));            % total transmitted power [W]
Nspan     = system.Nspan;
Loss      = system.Loss;
Hf_SNR    = system.mfil_snr;
Hf        = gpuArray(system.mfil);
set(t_sig,'POWER'     ,Pavg);
set(sig  ,'POWER'     ,Pavg);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                         Amplifier parameters                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ampli     = Ampliflat(Pavg,ch,ampli_par.G,ampli_par.e);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                         ESSFM PARAMETERS                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
options   = essfm.opt;
C0        = zeros(essfm.NC,1);
C0(1,1)   = 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                               TRAINING                                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    ux = gpuArray(complex(get(t_sig,'FIELDX')));
    uy = gpuArray(complex(get(t_sig,'FIELDY')));
    
    for i = 1:Nspan
        [ux, uy]      = ch.gpu_vectorial_ssfm(Pavg,t_sig,ux,uy);
        [ux, uy]      = ampli.gpu_AddNoise(t_sig,ux,uy);
    end
    
    fmin      = @(C) gpu_vec_essfm_opt(t_sig,ux,uy,dsp,C,Nspan,Loss,Hf);
    [C(1,:),err]=lsqnonlin(fmin,C0,[],[],options);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                      PROPAGATION                                        % 
%                    + BACKPROPAGATION ESSFM                              % 
%                    + SNR ESTIMATION                                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ux = gpuArray(complex(get(sig,'FIELDX')));
uy = gpuArray(complex(get(sig,'FIELDY')));
for i = 1:Nspan
    [ux, uy]      = ch.gpu_vectorial_ssfm(Pavg,sig,ux,uy);
    [ux, uy]      = ampli.gpu_AddNoise(sig,ux,uy);
end
ux_enh = ux;
uy_enh = uy;
sig_enh_rx = copy(sig);

if(dsp.nstep>01)
    for i = 1:Nspan
       [ux_enh, uy_enh]   = dsp.DBP_gpu_vec_essfm(Pavg*Loss,sig_enh_rx,C(1,:)',ux_enh,uy_enh);
    end
else
    for i = 1:round(Nspan*dsp.nstep)
        [ux_enh, uy_enh]  = dsp.DBP_gpu_vec_essfm(Pavg*Loss,sig_enh_rx,C(1,:)',ux_enh,uy_enh);
    end
end

set(sig_enh_rx,'FIELDX',gather(ux_enh));
set(sig_enh_rx,'FIELDY',gather(uy_enh));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           X-Pol                                    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

FIELDX_TX       = get(sig       ,'FIELDX_TX');
FIELDX_ENH_RX   = get(sig_enh_rx,'FIELDX'   );

FIELDX_ENH_RX   = ifft(fft(FIELDX_ENH_RX).*Hf_SNR);

rotx_enh        = angle(mean(FIELDX_ENH_RX.*conj(FIELDX_TX)));

FIELDX_ENH_RX   = FIELDX_ENH_RX*exp(-1i*rotx_enh);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           Y-Pol                                    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

FIELDY_TX       = get(sig       ,'FIELDY_TX');
FIELDY_ENH_RX   = get(sig_enh_rx,'FIELDY'   );

FIELDY_ENH_RX   = ifft(fft(FIELDY_ENH_RX).*Hf_SNR);

roty_enh        = angle(mean(FIELDY_ENH_RX.*conj(FIELDY_TX)));

FIELDY_ENH_RX   = FIELDY_ENH_RX*exp(-1i*roty_enh);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

snr_x = SNR(FIELDX_ENH_RX(1:sig.NT:end),FIELDX_TX(1:sig.NT:end));
snr_y = SNR(FIELDY_ENH_RX(1:sig.NT:end),FIELDY_TX(1:sig.NT:end));

avg_snr = 0.5*(snr_x+snr_y);

end

