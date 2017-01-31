classdef DSP
    %DSP Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        ch
        nstep;
        dz;
    end
    
    methods
        
        function obj = DSP(Channel,Ns_bprop)
            
            if isa(Channel,'Channel')
               obj.ch = Channel;
            else
                error(['Channel parameter has to be an instance',...
                    ' of Channel class']);
            end
            
            obj.nstep    = Ns_bprop;
            obj.dz       = Channel.Lf/Ns_bprop;

        end
        
        function [ux,uy]    = A2D (obj,Irx,sps)
            
            DecimationRate = sps / 2;                        %oversampling
            
            Irxdec = [  decimate( Irx(:,1), DecimationRate,16,'fir' ) ...
                        decimate( Irx(:,2), DecimationRate,16,'fir' ) ...
                        decimate( Irx(:,3), DecimationRate,16,'fir' ) ...
                        decimate( Irx(:,4), DecimationRate,16,'fir' ) ];
            
            ux = complex( Irxdec(:,1), Irxdec(:,2) );
            uy = complex( Irxdec(:,3), Irxdec(:,4) );
            
        end
        
        
        function sig        = DBP_scalar_ssfm (dsp,Pavg,sig)
            
            channel = dsp.ch;
            dz      = dsp.dz;
            
            omega = 2*pi*sig.SYMBOLRATE*sig.FN'*1e9;         % [rad/s]
            beta  = 0.5*omega.^2*channel.b2 + omega.^3*channel.b3/6;
            
            ux    = get(sig,'FIELDX');
                       
            if abs(channel.alphalin*dz) > 1e-6
                Leff     = -(1-exp(channel.alphalin*dz))/channel.alphalin;
            else
                Leff     = dz;
            end
            
            z    = dsp.dz*(0:dsp.nstep-1);
            
            xi   = -channel.gamma*Leff*exp(-channel.alphalin*z)*Pavg;
           
            halfdz      = dz/2;
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %                         HALF DZ GVD                        %
            %                              DZ SPM                        %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            ux = dsp.scalar_lin_step(-beta*halfdz,ux);
            ux = dsp.scalar_nl_step(xi,ux,1);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            for i=2:dsp.nstep  
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %                           DZ GVD                       %
                %                           DZ SPM                       %
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                ux = dsp.scalar_lin_step(-beta*dz,ux);
                ux = dsp.scalar_nl_step(xi,ux,i);
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
               
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %                      LAST HALF DZ GVD                      %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            ux = dsp.scalar_lin_step(-beta*halfdz,ux);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            set(sig,'FIELDX',ux);
            
        end
        
        function sig        = DBP_vec_ssfm (dsp,Pavg,sig)
            channel = dsp.ch;
            dz      = dsp.dz;
            
            omega = 2*pi*sig.SYMBOLRATE*sig.FN'*1e9;         % [rad/s]
            beta  = 0.5*omega.^2*channel.b2 + omega.^3*channel.b3/6;
            
            ux    = get(sig,'FIELDX');
            uy    = get(sig,'FIELDY');
                        
            if abs(channel.alphalin*dz) > 1e-6
                Leff     = -(1-exp(channel.alphalin*dz))/channel.alphalin;
            else
                Leff     = dz;
            end
            
            z    = dsp.dz*(0:dsp.nstep-1);
            
            xi   = -channel.gamma*Leff*exp(-channel.alphalin*z)*Pavg;
           
            halfdz      = dz/2;
                
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %                           HALF DZ GVD                      %
            %                                DZ SPM                      %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            [ux,uy] = dsp.vec_lin_step(-beta*halfdz,ux,uy);
            [ux,uy] = dsp.vec_nl_step(xi(1),ux,uy);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        
            for i=2:dsp.nstep
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %                           DZ GVD                       %
                %                           DZ SPM                       %
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                [ux,uy] = dsp.vec_lin_step(-beta*dz,ux,uy); 
                [ux,uy] = dsp.vec_nl_step(xi(i),ux,uy); 
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
               
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %                      LAST HALF DZ GVD                      %
            %                      LAST      DZ SPM                      %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            [ux,uy] = dsp.vec_lin_step(-beta*halfdz,ux,uy);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            set(sig,'FIELDX',ux);
            set(sig,'FIELDY',uy);
        end
        
        function sig        = DBP_essfm(dsp,Pavg,sig,C)

            channel = dsp.ch;
            dz = dsp.dz;
            
            omega = 2*pi*sig.SYMBOLRATE*sig.FN'*1e9;         % [rad/s]
            beta  = 0.5*omega.^2*channel.b2 + omega.^3*channel.b3/6;
     
            ux    = get(sig,'FIELDX');
            
            if abs(channel.alphalin*dz) > 1e-6
                Leff     = -(1-exp(channel.alphalin*dz))/channel.alphalin;
            else
                Leff     = dz;
            end
            
            z    = dsp.dz*(0:dsp.nstep-1);
            
            xi   = -channel.gamma*Leff*exp(-channel.alphalin*z)*Pavg;
           
            halfdz      = dz/2;
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %                         HALF DZ GVD                        %
            %                              DZ SPM                        %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            ux = dsp.scalar_lin_step(-beta*halfdz,ux);
            ux = dsp.scalar_essfm_nl_step(-xi(1)*C,ux);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            for i=2:dsp.nstep  
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %                           DZ GVD                       %
                %                           DZ SPM                       %
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                ux = dsp.scalar_lin_step(-beta*dz,ux);
                ux = dsp.scalar_essfm_nl_step(-xi(i)*C,ux);
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
               
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %                      LAST HALF DZ GVD                      %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            ux = dsp.scalar_lin_step(-beta*halfdz,ux);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            set(sig,'FIELDX',ux);
        end
        
        function sig        = DBP_vec_essfm (dsp,Pavg,sig,C)
            
            channel = dsp.ch;
            dz = dsp.dz;
            
            omega = 2*pi*sig.SYMBOLRATE*sig.FN'*1e9;         % [rad/s]
            beta  = 0.5*omega.^2*channel.b2 + omega.^3*channel.b3/6;
     
            ux    = get(sig,'FIELDX');
            uy    = get(sig,'FIELDY');
            
            if abs(channel.alphalin*dz) > 1e-6
                Leff     = -(1-exp(channel.alphalin*dz))/channel.alphalin;
            else
                Leff     = dz;
            end
            
            z    = dsp.dz*(0:dsp.nstep-1);
            
            xi   = -channel.gamma*Leff*exp(-channel.alphalin*z)*Pavg;
           
            halfdz      = dz/2;
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %                           HALF DZ GVD                      %
            %                                DZ SPM                      %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            [ux,uy] = dsp.vec_lin_step(-beta*halfdz,ux,uy);;
            [ux,uy] = dsp.vec_nl_essfm_step(-xi(1)*C,ux,uy);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        
            for i=2:dsp.nstep  
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %                           DZ GVD                       %
                %                           DZ SPM                       %
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                [ux,uy] = dsp.vec_lin_step(-beta*dz,ux,uy);
                [ux,uy] = dsp.vec_nl_essfm_step(-xi(i)*C,ux,uy);
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %                      LAST HALF DZ GVD                      %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            [ux,uy] = dsp.vec_lin_step(-beta*halfdz,ux,uy);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            set(sig,'FIELDX',ux);
            set(sig,'FIELDY',uy);
        end
                       
    end
    
    methods (Access = private)
                
        function ux         = scalar_lin_step(obj,betaxdz,ux)
            
            Hf = exp(-1i*betaxdz);
            ux = ifft( fft(ux) .* Hf);
            
        end
        
        function ux         = scalar_nl_step(obj,xi,ux,step)
            
            pow = real(ux).^2 + imag(ux).^2;
            ux = ux.*exp(-1i*xi(step).*pow);
            
        end
        
        function [ux,uy]    = vec_lin_step(obj,betaxdz,ux,uy)
            
            Hf = exp(-1i*betaxdz);
            
            ux = ifft( fft(ux) .* Hf);
            uy = ifft( fft(uy) .* Hf);
            
        end
        
        function [ux,uy]    = vec_nl_step(obj,xi,ux,uy)
            
            pow = real(ux).^2 + imag(ux).^2 + real(uy).^2 + imag(uy).^2;
            ux = ux.*exp(-1i*xi.*pow);
            uy = uy.*exp(-1i*xi.*pow);
            
        end
        
        function [ux]       = scalar_essfm_nl_step(obj,C,ux)
            
            pow = real(ux).^2 + imag(ux).^2 ;
            
            M=length(pow);
            N=length(C);
            per_pow =[pow(M-N+2:M);pow;pow(1:N-1)];
            
            CC=[flipud(C);C(2:end)];
            theta=conv(per_pow,CC,'valid');
            
            ux = ux.*exp(1i*theta);
            
        end
        
        function [ux,uy]    = vec_nl_essfm_step(obj,C,ux,uy)
            
            pow = real(ux).^2 + imag(ux).^2 + real(uy).^2 + imag(uy).^2;
            
            M=length(pow);
            N=length(C);
            per_pow =[pow(M-N+2:M);pow;pow(1:N-1)];
            
            CC=[flipud(C);C(2:end)];
            theta=conv(per_pow,CC,'valid');
            
            ux = ux.*exp(1i*theta);
            uy = uy.*exp(1i*theta);
            
        end
        
    end
    
end
