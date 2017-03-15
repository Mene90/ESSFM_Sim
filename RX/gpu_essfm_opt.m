function [ f ] = gpu_essfm_opt( sig,ux,dsp,C,Nspan,Loss,Hf )
%GPU_ESSFM_OPT Summary of this function goes here
%   Detailed explanation goes here
        rx_sig     = copy(sig);
              
        NT      = get(rx_sig,'NT'       );
        Pavg    = get(rx_sig,'POWER'    );
        
        if(dsp.nstep>=1)
            for i = 1:Nspan
                ux = dsp.DBP_gpu_essfm(Pavg*Loss,rx_sig,C,ux);
            end
        else
            for i=1:round(Nspan*dsp.nstep)
                ux = dsp.DBP_gpu_essfm(Pavg*Loss,rx_sig,C,ux);
            end
        end
        
        ux_out  = get(rx_sig,'FIELDX_TX');   
        
        ux    = ifft(fft(ux).*Hf);
        
        ux_rx   = gather(ux);
        
        ux_rx   = ux_rx(1:NT:end);
%         ux_out  = ux_out(1:NT:end);
        
        rot=angle(mean(ux_rx.*conj(ux_out)));
        ux_rx=ux_rx*exp(-1i*rot);
        
        f = [real(ux_rx-ux_out).'; imag(ux_rx-ux_out).'];
        
        return

end

