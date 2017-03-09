function [ f ] = gpu_vec_essfm_opt( sig,ux,uy,dsp,C,Nspan,Loss,Hf )
%VEC_ESSFM_OPT Summary of this function goes here
%   Detailed explanation goes here
        rx_sig     = copy(sig);
              
        NT      = get(rx_sig,'NT'       );
        Pavg    = get(rx_sig,'POWER'    );
        
        if(dsp.nstep>=1)
            for i = 1:Nspan
                [ux, uy] = dsp.DBP_gpu_vec_essfm(Pavg*Loss,rx_sig,C,ux,uy);
            end
        else
            for i=1:round(Nspan*dsp.nstep)
                [ux, uy] = dsp.DBP_gpu_vec_essfm(Pavg*Loss,rx_sig,C,ux,uy);
            end
        end
        
        
        ux_out  = get(rx_sig,'FIELDX_TX');
        uy_out  = get(rx_sig,'FIELDY_TX');  

        
        ux    = ifft(fft(ux).*Hf);
        uy    = ifft(fft(uy).*Hf);
        
        ux_rx   = gather(ux);
        uy_rx   = gather(uy);
        
        ux_rx   = ux_rx(1:NT:end);
        uy_rx   = uy_rx(1:NT:end);
        ux_out  = ux_out(1:NT:end);
        uy_out  = uy_out(1:NT:end);
        
        rotx    = angle(mean(ux_rx.*conj(ux_out)));
        roty    = angle(mean(uy_rx.*conj(uy_out)));
        
        ux_rx   = ux_rx*exp(-1i*rotx);
        uy_rx   = uy_rx*exp(-1i*roty);
        
        f = [real(ux_rx-ux_out); imag(ux_rx-ux_out);...
             real(uy_rx-uy_out); imag(uy_rx-uy_out);];
   
        return


end

