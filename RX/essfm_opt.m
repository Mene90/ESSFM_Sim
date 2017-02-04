function  f  = essfm_opt( sig,dsp,C,Nspan,Loss, Hf )
%ESSFM_OPT Summary of this function goes here
%   Detailed explanation goes here
        
        rx_sig     = copy(sig);
              
        NT      = get(rx_sig,'NT'       );
        Pavg    = get(rx_sig,'POWER'    );
        
        for i = 1:Nspan
            rx_sig = dsp.DBP_essfm(Pavg*Loss,rx_sig,C);
        end
        
        ux_out  = get(rx_sig,'FIELDX_TX');  
        ux_rx   = get(rx_sig,'FIELDX'   );
        
        
        ux_rx    = ifft(fft(ux_rx).*Hf);
        
        ux_rx   = ux_rx(1:NT:end);
        ux_out  = ux_out(1:NT:end);
        
        rot=angle(mean(ux_rx.*conj(ux_out)));
        ux_rx=ux_rx*exp(-1i*rot);
        
        f = [real(ux_rx-ux_out); imag(ux_rx-ux_out);];
        
        return
end

