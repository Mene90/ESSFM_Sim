classdef DSP
    %DSP Summary of this class goes here
    %    Detailed explanation goes here
    
    properties
        ch
        nstep;
        dz;
    end
    
    methods
        
        function obj        = DSP(Channel,Ns_bprop)
            
            if isa(Channel,'Channel')
               obj.ch = Channel;
            else
                error(['Channel parameter has to be an instance',...
                    ' of Channel class']);
            end
            
            obj.nstep    = Ns_bprop;
            obj.dz       = Channel.Lf/Ns_bprop;

        end
        
        function ux         = DBP_gpu_scalar_ssfm (dsp,Pavg,sig,ux)
            
            channel = dsp.ch;
            dz      = dsp.dz;
            halfdz   = dz/2;
            
            omega = 2*pi*sig.SYMBOLRATE*sig.FN'*1e9;         % [rad/s]
            betaz  = gpuArray(complex((0.5*omega.^2*channel.b2 + omega.^3*channel.b3/6)*dz));
            betahz  = gpuArray(complex((0.5*omega.^2*channel.b2 + omega.^3*channel.b3/6)*halfdz));
               
            if(dsp.nstep>=1)
                
                if abs(-channel.alphalin*dz) > 1e-6
                    Leff     = -(1-exp(channel.alphalin*dz))/channel.alphalin;
                else
                    Leff     = dz;
                end
                
                z    = dsp.dz*(0:dsp.nstep-1);
                
                xi   = gpuArray(-channel.gamma*Leff*exp(channel.alphalin*z)*Pavg);
                
                                
            else
                
                
                Leff     = (1.0-exp(channel.alphalin*channel.Lf))/(-channel.alphalin*channel.Lf);
                xi       = gpuArray(-channel.gamma*Leff*Pavg*dz);
                
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %                         HALF DZ GVD                        %
            %                              DZ SPM                        %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            ux = dsp.scalar_lin_step(-betahz,ux);
            ux = dsp.scalar_nl_step(xi(1),ux);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            for i=2:dsp.nstep  
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %                           DZ GVD                       %
                %                           DZ SPM                       %
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                ux = dsp.scalar_lin_step(-betaz,ux);
                ux = dsp.scalar_nl_step(xi(i),ux);
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
               
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %                      LAST HALF DZ GVD                      %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            ux = dsp.scalar_lin_step(-betahz,ux);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        
        end
        
        function [ux, uy]   = DBP_gpu_vec_ssfm (dsp,Pavg,sig,ux,uy)
                            
            channel = dsp.ch;
            dz      = dsp.dz;
            halfdz      = dz/2;
            
            omega = 2*pi*sig.SYMBOLRATE*sig.FN'*1e9;         % [rad/s]
            betaz  = gpuArray(complex((0.5*omega.^2*channel.b2 + omega.^3*channel.b3/6)*dz));
            betahz  = gpuArray(complex((0.5*omega.^2*channel.b2 + omega.^3*channel.b3/6)*halfdz));
            
            
            if(dsp.nstep>=1)
                        
                if abs(channel.alphalin*dz) > 1e-6
                    Leff     = -(1-exp(channel.alphalin*dz))/channel.alphalin;
                else
                    Leff     = dz;
                end
                
                z    = dsp.dz*(0:dsp.nstep-1);
                
                xi   = gpuArray(-channel.gamma*Leff*exp(channel.alphalin*z)*Pavg);
                                
            else
                
                Leff     = (1.0-exp(channel.alphalin*channel.Lf))/(-channel.alphalin*channel.Lf);
                xi       = gpuArray(-channel.gamma*Leff*Pavg*dz);
            
            end
                
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %                           HALF DZ GVD                      %
            %                                DZ SPM                      %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            [ux,uy] = dsp.vec_lin_step(-betahz,ux,uy);
            [ux,uy] = dsp.vec_nl_step(xi(1),ux,uy);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        
            for i=2:dsp.nstep
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %                           DZ GVD                       %
                %                           DZ SPM                       %
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                [ux,uy] = dsp.vec_lin_step(-betaz,ux,uy); 
                [ux,uy] = dsp.vec_nl_step(xi(i),ux,uy); 
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
               
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %                      LAST HALF DZ GVD                      %
            %                      LAST      DZ SPM                      %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            [ux,uy] = dsp.vec_lin_step(-betahz,ux,uy);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
        end
        
        function [ux, uy]   = DBP_gpu_vec_ssfm_disp_comp(dsp,Pavg,sig,ux,uy)
                            
            channel = dsp.ch;
            dz      = dsp.dz;
            halfdz      = dz/2;
            
            omega = 2*pi*sig.SYMBOLRATE*sig.FN'*1e9;         % [rad/s]
            betaz  = gpuArray(complex((0.5*omega.^2*channel.b2 + omega.^3*channel.b3/6)*dz));
            betahz  = gpuArray(complex((0.5*omega.^2*channel.b2 + omega.^3*channel.b3/6)*halfdz));
            
            
            if(dsp.nstep>=1)
                        
                if abs(channel.alphalin*dz) > 1e-6
                    Leff     = -(1-exp(channel.alphalin*dz))/channel.alphalin;
                else
                    Leff     = dz;
                end
                
                z    = dsp.dz*(0:dsp.nstep-1);
                
                xi   = gpuArray(-channel.gamma*Leff*exp(channel.alphalin*z)*Pavg);
                                
            else
                
                Leff     = (1.0-exp(channel.alphalin*channel.Lf))/(-channel.alphalin*channel.Lf);
                xi       = gpuArray(-channel.gamma*Leff*Pavg*dz);
            
            end
                
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %                           HALF DZ GVD                      %
            %                                DZ SPM                      %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            [ux,uy] = dsp.vec_lin_step(-betahz,ux,uy);
%             [ux,uy] = dsp.vec_nl_step(xi(1),ux,uy);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        
            for i=2:dsp.nstep
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %                           DZ GVD                       %
                %                           DZ SPM                       %
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                [ux,uy] = dsp.vec_lin_step(-betaz,ux,uy); 
%                 [ux,uy] = dsp.vec_nl_step(xi(i),ux,uy); 
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
               
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %                      LAST HALF DZ GVD                      %
            %                      LAST      DZ SPM                      %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            [ux,uy] = dsp.vec_lin_step(-betahz,ux,uy);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
        end
        
        function ux         = DBP_gpu_essfm(dsp,Pavg,sig,C,ux)
            
            channel = dsp.ch;
            dz = dsp.dz;
            halfdz      = dz/2;
            
            omega = 2*pi*sig.SYMBOLRATE*sig.FN'*1e9;         % [rad/s]
            betaz  = gpuArray(complex((0.5*omega.^2*channel.b2 + omega.^3*channel.b3/6)*dz));
            betahz = gpuArray(complex((0.5*omega.^2*channel.b2 + omega.^3*channel.b3/6)*halfdz));
            
                        
            if(dsp.nstep>=1)
                if abs(channel.alphalin*dz) > 1e-6
                    Leff     = -(1-exp(channel.alphalin*dz))/channel.alphalin;
                else
                    Leff     = dz;
                end
                
                z    = dsp.dz*(0:dsp.nstep-1);
                
                xi   = gpuArray(-channel.gamma*Leff*exp(channel.alphalin*z)*Pavg); 
                
            else
                
                Leff     = (1.0-exp(channel.alphalin*channel.Lf))/(-channel.alphalin*channel.Lf);
                xi       = gpuArray(-channel.gamma*Leff*Pavg*dz);
                
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %                         HALF DZ GVD                        %
            %                              DZ SPM                        %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            ux = dsp.scalar_lin_step(-betahz,ux);
            ux = dsp.scalar_essfm_nl_step(-xi(1)*C,ux);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            for i=2:dsp.nstep
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %                           DZ GVD                       %
                %                           DZ SPM                       %
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                ux = dsp.scalar_lin_step(-betaz,ux);
                ux = dsp.scalar_essfm_nl_step(-xi(i)*C,ux);
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %                      LAST HALF DZ GVD                      %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            ux = dsp.scalar_lin_step(-betahz,ux);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
        end
        
        function [ux,uy]    = DBP_gpu_vec_essfm (dsp,Pavg,sig,C,ux,uy)
            
            channel = dsp.ch;
            dz = dsp.dz;
            halfdz      = dz/2;
            
            omega = 2*pi*sig.SYMBOLRATE*sig.FN'*1e9;         % [rad/s]
            betaz  = gpuArray(complex((0.5*omega.^2*channel.b2 + omega.^3*channel.b3/6)*dz));
            betahz  = gpuArray(complex((0.5*omega.^2*channel.b2 + omega.^3*channel.b3/6)*halfdz));
                        
            if(dsp.nstep>1)
                
                if abs(channel.alphalin*dz) > 1e-6
                    Leff     = -(1-exp(channel.alphalin*dz))/channel.alphalin;
                else
                    Leff     = dz;
                end
                
                z    = dsp.dz*(0:dsp.nstep-1);
                
                xi   = gpuArray(-channel.gamma*Leff*exp(channel.alphalin*z)*Pavg);
                                                
            else
               
                Leff     = (1.0-exp(channel.alphalin*channel.Lf))/(-channel.alphalin*channel.Lf)*dz;
                xi       = gpuArray(-channel.gamma*Leff*Pavg);
                
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %                           HALF DZ GVD                      %
            %                                DZ SPM                      %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            [ux,uy] = dsp.vec_lin_step(-betahz,ux,uy);
            [ux,uy] = dsp.vec_nl_essfm_step(-xi(1)*C,ux,uy);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            for i=2:dsp.nstep
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %                           DZ GVD                       %
                %                           DZ SPM                       %
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                [ux,uy] = dsp.vec_lin_step(-betaz,ux,uy);
                [ux,uy] = dsp.vec_nl_essfm_step(-xi(i)*C,ux,uy);
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %                      LAST HALF DZ GVD                      %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            [ux,uy] = dsp.vec_lin_step(-betahz,ux,uy);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
        end
       
        function sig        = DBP_vec_ssfm_disp_comp(dsp,Pavg,sig)
        channel = dsp.ch;
            dz      = dsp.dz;
            
            omega = 2*pi*sig.SYMBOLRATE*sig.FN'*1e9;         % [rad/s]
            beta  = 0.5*omega.^2*channel.b2 + omega.^3*channel.b3/6;
            
            ux    = get(sig,'FIELDX');
            uy    = get(sig,'FIELDY');
            if(dsp.nstep>=1)
                
                if abs(channel.alphalin*dz) > 1e-6
                    Leff     = -(1-exp(channel.alphalin*dz))/channel.alphalin;
                else
                    Leff     = dz;
                end
                
                z    = dsp.dz*(0:dsp.nstep-1);
                
                xi   = -channel.gamma*Leff*exp(channel.alphalin*z)*Pavg;
                
                halfdz      = dz/2;
            else
                
                halfdz   = dz/2;
                Leff     = (1.0-exp(channel.alphalin*channel.Lf))/(-channel.alphalin*channel.Lf);
                xi       = -channel.gamma*Leff*Pavg*dz;
                
            end
                
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %                           HALF DZ GVD                      %
            %                                DZ SPM                      %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            [ux,uy] = dsp.vec_lin_step(-beta*halfdz,ux,uy);
%             [ux,uy] = dsp.vec_nl_step(xi(1),ux,uy);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        
            for i=2:dsp.nstep
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %                           DZ GVD                       %
                %                           DZ SPM                       %
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                [ux,uy] = dsp.vec_lin_step(-beta*dz,ux,uy); 
%                 [ux,uy] = dsp.vec_nl_step(xi(i),ux,uy); 
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
                
        function sig        = DBP_scalar_ssfm (dsp,Pavg,sig)
            
            channel = dsp.ch;
            dz      = dsp.dz;
            
            omega = 2*pi*sig.SYMBOLRATE*sig.FN'*1e9;         % [rad/s]
            beta  = 0.5*omega.^2*channel.b2 + omega.^3*channel.b3/6;
            
            ux    = get(sig,'FIELDX');
               
            if(dsp.nstep>=1)
                
                if abs(-channel.alphalin*dz) > 1e-6
                    Leff     = -(1-exp(channel.alphalin*dz))/channel.alphalin;
                else
                    Leff     = dz;
                end
                
                z    = dsp.dz*(0:dsp.nstep-1);
                
                xi   = -channel.gamma*Leff*exp(channel.alphalin*z)*Pavg;
                
                halfdz      = dz/2;
                
            else
                
                halfdz   = dz/2;
                Leff     = (1.0-exp(channel.alphalin*channel.Lf))/(-channel.alphalin*channel.Lf);
                xi       = -channel.gamma*Leff*Pavg*dz;
                
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %                         HALF DZ GVD                        %
            %                              DZ SPM                        %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            ux = dsp.scalar_lin_step(-beta*halfdz,ux);
            ux = dsp.scalar_nl_step(xi(1),ux);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            for i=2:dsp.nstep  
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %                           DZ GVD                       %
                %                           DZ SPM                       %
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                ux = dsp.scalar_lin_step(-beta*dz,ux);
                ux = dsp.scalar_nl_step(xi(i),ux);
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
            
            if(dsp.nstep>=1)
                        
                if abs(channel.alphalin*dz) > 1e-6
                    Leff     = -(1-exp(channel.alphalin*dz))/channel.alphalin;
                else
                    Leff     = dz;
                end
                
                z    = dsp.dz*(0:dsp.nstep-1);
                
                xi   = -channel.gamma*Leff*exp(channel.alphalin*z)*Pavg;
                
                halfdz      = dz/2;
                
            else
                
                halfdz   = dz/2;
                Leff     = (1.0-exp(channel.alphalin*channel.Lf))/(-channel.alphalin*channel.Lf);
                xi       = -channel.gamma*Leff*Pavg*dz;
            
            end
                
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
            
            if(dsp.nstep>=1)
                if abs(channel.alphalin*dz) > 1e-6
                    Leff     = -(1-exp(channel.alphalin*dz))/channel.alphalin;
                else
                    Leff     = dz;
                end
                
                z    = dsp.dz*(0:dsp.nstep-1);
                
                xi   = -channel.gamma*Leff*exp(channel.alphalin*z)*Pavg;
                
                halfdz      = dz/2;
            else
                
                halfdz   = dz/2;
                Leff     = (1.0-exp(channel.alphalin*channel.Lf))/(-channel.alphalin*channel.Lf);
                xi       = -channel.gamma*Leff*Pavg*dz;
                
            end
            
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
            
            if(dsp.nstep>=1)
                
                if abs(channel.alphalin*dz) > 1e-6
                    Leff     = -(1-exp(channel.alphalin*dz))/channel.alphalin;
                else
                    Leff     = dz;
                end
                
                z    = dsp.dz*(0:dsp.nstep-1);
                
                xi   = -channel.gamma*Leff*exp(channel.alphalin*z)*Pavg;
                
                halfdz      = dz/2;
                
            else
                
                halfdz   = dz/2;
                Leff     = (1.0-exp(channel.alphalin*channel.Lf))/(-channel.alphalin*channel.Lf);
                xi       = -channel.gamma*Leff*Pavg*dz;
                
            end
            
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
        
        function ux         = scalar_nl_step(obj,xi,ux)
            
            pow = abs(ux).^2;
            ux = ux.*exp(-1i*xi*pow);
            
        end
        
        function [ux,uy]    = vec_lin_step(obj,betaxdz,ux,uy)
            
            Hf = exp(-1i*betaxdz);
            
            ux = ifft( fft(ux) .* Hf);
            uy = ifft( fft(uy) .* Hf);
            
        end
        
        function [ux,uy]    = vec_nl_step(obj,xi,ux,uy)
            
            pow = abs(ux).^2+abs(uy).^2;
            ux = ux.*exp(-1i*xi.*pow);
            uy = uy.*exp(-1i*xi.*pow);
            
        end
        
        function [ux]       = scalar_essfm_nl_step(obj,C,ux)
            
            
            pow = real(ux).^2 + imag(ux).^2 ;
            
            M=length(pow);
            N=length(C);
            per_pow =[pow(M-N+2:M);pow;pow(1:N-1)];
            if(length(C)>1)
                CC=[flipud(C);C(2:end)];
            else
                CC = [flipud(C)];
            end
            theta=conv(per_pow,CC,'valid');
            
            ux = ux.*exp(1i*theta);
            
        end
        
        function [ux,uy]    = vec_nl_essfm_step(obj,C,ux,uy)
            
            
            pow = abs(ux).^2 + abs(uy).^2;
            
            M=length(pow);
            N=length(C);
            per_pow =[pow(M-N+2:M);pow;pow(1:N-1)];
            if(length(C)>1)
                CC = [flipud(C);C(2:end)];
            else
                CC = [flipud(C)];
            end
            theta=conv(per_pow,CC,'valid');
            
            ux = ux.*exp(1i*theta*8/9);
            uy = uy.*exp(1i*theta*8/9);
            
        end
        
    end
    
end

































%             M=length(ux);
%             N=length(C);
%             
%             % Alloca spazio per il segnale
% %             y=zeros(M,1);
%             % Calcola moduli quadri (con prolungamento periodico)
%             xx=abs(ux).^2 + abs(uy).^2;
%             xx=[xx(M-N+2:M);xx;xx(1:N-1)];    %periodicamente
%             %xx=[0;xx;0];                     %con zeri
%             
%             % Duplica il vettore dei coefficienti per avere risposta simmetrica
%             CC=[flipud(C);C(2:end)];
%             
%             % Filtra il vettore dei moduli quadri con i coefficienti del filtro dato:
%             theta=conv(xx,CC,'valid');
%             
%             % Calcola il segnale in uscita:
%             ux=ux.*exp(1i*theta);
%             uy=uy.*exp(1i*theta);


%             M=length(ux);
%             N=length(C);
%             
%             % Alloca spazio per il segnale
% %             y=zeros(M,1);
%             % Calcola moduli quadri (con prolungamento periodico)
%             xx=abs(ux).^2;
%             xx=[xx(M-N+2:M);xx;xx(1:N-1)];    %periodicamente
%             %xx=[0;xx;0];                     %con zeri
%             
%             % Duplica il vettore dei coefficienti per avere risposta simmetrica
%             CC=[flipud(C);C(2:end)];
%             
%             % Filtra il vettore dei moduli quadri con i coefficienti del filtro dato:
%             theta=conv(xx,CC,'valid');
%             
%             % Calcola il segnale in uscita:
%             ux=ux.*exp(1i*theta);

