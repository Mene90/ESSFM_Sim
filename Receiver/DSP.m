classdef DSP < handle & matlab.mixin.SetGet
    %DSP Summary of this class goes here
    %    Detailed explanation goes here
    
    properties  
        ch;
        nstep;
        dz;
        C;
    end
    
    methods     
        
        function dsp  = DSP(Channel,Ns_bprop)       
            
            if isa(Channel,'Channel')
                dsp.ch = Channel;
            else
                error(['Channel parameter has to be an instance',...
                    ' of Channel class']);
            end
            
            dsp.nstep = Ns_bprop;
            dsp.dz    = Channel.Lf/Ns_bprop;

        end
        
        function setessfmcoeff(dsp,NC,Pavg,Nsymb,Nt,symbrate,pls,ampli,Nspan)
            
            Loss    = 10^(-(dsp.ch.alphadB*dsp.ch.Lf*1e-3)*0.1);
            C0      = zeros(NC,1);
            C0(1,1) = 1;
            options = optimset('Algorithm','trust-region-reflective',...
                               'Display','off','Jacobian','off',...
                               'DerivativeCheck','off','TolFun',...
                                1e-13,'TolX',1e-13);
                   
                            
            sig = Signal(Nsymb,Nt,symbrate);
            Hf  = transpose(filt(pls,sig.FN));
            
           [patx(:,1), patmatx] = Pattern.debruijn(1,4,Nsymb);
           [paty(:,1), patmaty] = Pattern.debruijn(2,4,Nsymb);
           
           E                    = Laser.GetLaserSource(Pavg, Nsymb * Nt);
           
           set(sig,'POWER'    ,Pavg);
           set(sig,'FIELDX'   ,Modulator.ApplyModulation(E, 2*patmatx-1, sig, pls));
           set(sig,'FIELDY'   ,Modulator.ApplyModulation(E, 2*patmaty-1, sig, pls));
           set(sig,'FIELDX_TX',Modulator.ApplyModulation(E, 2*patmatx-1, sig, pls));
           set(sig,'FIELDY_TX',Modulator.ApplyModulation(E, 2*patmaty-1, sig, pls));
           
           
           dsp.ch.propagate(Nspan,ampli,sig);
           
           fmin = @(C) essfm_coeff_opt(sig,dsp,C,Nspan,Loss,Hf);
           [dsp.C,err] = lsqnonlin(fmin,C0,[],[],options);
        end
        
        function backpropagation(dsp,Pavg,sig,Nspan,type,comp_ch)
            r     = 0.5;           
            
            omega = 2*pi*sig.SYMBOLRATE*sig.FN'*1e9;             % [rad/s]
            betaz = gpuArray(complex((0.5 * omega.^2*dsp.ch.b2...
                                          + omega.^3*dsp.ch.b3/6)*dsp.dz));

            dx    = (dsp.ch.Lf>=dsp.dz)*dsp.dz     +... 
                    (1-(dsp.ch.Lf>=dsp.dz))*dsp.ch.Lf;
                
                
            k     = (dsp.nstep<1)*dsp.dz/dsp.ch.Lf +...
                    (1-(dsp.nstep<1))*1;
            
            Leff  = dsp.dz;
            if abs(dsp.ch.alphalin*dsp.dz) > 1e-6
                Leff = -k*(1-exp(dsp.ch.alphalin*dx))/dsp.ch.alphalin;
            end
            
            z     = dsp.dz *(0:ceil(dsp.nstep)-1);
            xi    = gpuArray(-dsp.ch.gamma*Leff*exp(dsp.ch.alphalin*z)*Pavg);
            
            if(any(sig.FIELDY))
                xi = xi*8/9;
            end

            fiberpartion = (dsp.nstep<1)*1+ (1-(dsp.nstep<1))*Nspan;
            
            set(sig,'FIELDX', gpuArray(complex(get(sig,'FIELDX'))));
            set(sig,'FIELDY', gpuArray(complex(get(sig,'FIELDY'))));
            
            if strcmpi(type,'ssfm')
                for i=1:fiberpartion
                    ssfm(dsp,sig,xi,Nspan,betaz,r);
                end
            elseif strcmpi(type,'essfm')
                for i=1:fiberpartion
                    essfm(dsp,sig,xi,Nspan,betaz,r);
                end
            else
                error('backpropagation method not recognized');
            end
            
            set(sig,'FIELDX', gather(get(sig,'FIELDX')));
            set(sig,'FIELDY', gather(get(sig,'FIELDY')));
        
        end
                
        function matchedfilter(~,sig,Hf)            
            set(sig,'FIELDX', ifft(fft(sig.FIELDX).*Hf));
            set(sig,'FIELDY', ifft(fft(sig.FIELDY).*Hf));
        end
        
        function downsampling(~,sig)
            if(size(sig.FIELDX,1) == size(sig.FIELDX_TX,1)*sig.NT)
                set(sig,'FIELDX',sig.FIELDX(1:sig.NT:end));
                set(sig,'FIELDY',sig.FIELDY(1:sig.NT:end));
            end
        end        
        
        function nlpnmitigation(~,sig) 
            
            rotx = angle(mean(sig.FIELDX.*conj(get(sig,'FIELDX_TX'))));
            roty = angle(mean(sig.FIELDY.*conj(get(sig,'FIELDY_TX'))));
            
            set(sig,'FIELDX', sig.FIELDX*exp(-1i*rotx));
            set(sig,'FIELDY', sig.FIELDY*exp(-1i*roty));
            
        end
        
    end
    
    methods (Access = private)   
       
        function ssfm(dsp,sig,xi,Nspan,betaz,r)    
            steps = dsp.nstep;
            if(dsp.nstep < 1)
                steps = Nspan * dsp.ch.Lf/dsp.dz;
                xi    = ones(1,steps)*xi;
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %                           HALF DZ GVD                      %
            %                                DZ SPM                      %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            lin_step(dsp,-betaz*(1-r),sig);
            nl_ssfm_step(dsp,-xi(1),sig);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            for i=2:steps
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %                           DZ GVD                       %
                %                           DZ SPM                       %
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                lin_step(dsp,-betaz,sig);
                nl_ssfm_step(dsp,-xi(i),sig);
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %                      LAST HALF DZ GVD                      %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            lin_step(dsp,-betaz*r,sig);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        end
        
        function essfm(dsp,sig,xi,Nspan,betaz,r)   
             
            steps = dsp.nstep;
            if(dsp.nstep <= 1)
                steps = Nspan * dsp.ch.Lf/dsp.dz;
                xi    = ones(1,steps)*xi;
            end          
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %                           HALF DZ GVD                      %
            %                                DZ SPM                      %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            lin_step(dsp,-betaz*(1-r),sig);
            nl_essfm_step(dsp,-xi(1)*dsp.C,sig);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%            
            for i=2:steps
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %                           DZ GVD                       %
                %                           DZ SPM                       %
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                lin_step(dsp,-betaz,sig);
                nl_essfm_step(dsp,-xi(i)*dsp.C,sig);
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            end            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %                      LAST HALF DZ GVD                      %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            lin_step(dsp,-betaz*r,sig);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%            
        end
        
        function nl_ssfm_step(~,xi,sig)            
            pow = abs(sig.FIELDX).^2 + abs(sig.FIELDY).^2;
            set(sig,'FIELDX', sig.FIELDX.*exp(1i*xi.*pow));            
            set(sig,'FIELDY', sig.FIELDY.*exp(1i*xi.*pow));             
        end
        
        function nl_essfm_step(~,C,sig)            

            pow = abs(sig.FIELDX).^2 + abs(sig.FIELDY).^2;
            M   = length(pow);
            N   = length(C);
            per_pow =[pow(M-N+2:M);pow;pow(1:N-1)];
            
            if(length(C)>1)
                CC = [flipud(C);C(2:end)];
            else
                CC = [flipud(C)];
            end
            
            theta = conv(per_pow,CC,'valid');
            set(sig,'FIELDX', sig.FIELDX.*exp(1i*theta*8/9));
            set(sig,'FIELDY', sig.FIELDY.*exp(1i*theta*8/9)); 
            
        end
        
        function lin_step(~,betaxdz,sig)           
            set(sig,'FIELDX',ifft( fft(sig.FIELDX).* exp(-1i*betaxdz)));
            set(sig,'FIELDY',ifft( fft(sig.FIELDY).* exp(-1i*betaxdz)));             
        end
        
    end
end

function [ f ] = essfm_coeff_opt(sig,dsp,C,Nspan,Loss,Hf)  

        rx_sig = copy(sig);
        set(dsp,'C',C);
        dsp.backpropagation(get(rx_sig,'POWER')*Loss,rx_sig,Nspan,'essfm');
        dsp.matchedfilter(rx_sig,Hf);
        dsp.nlpnmitigation(rx_sig);
        
        X = [real(get(rx_sig,'FIELDX')-get(rx_sig,'FIELDX_TX')).';...
             imag(get(rx_sig,'FIELDX')-get(rx_sig,'FIELDX_TX')).'];
        Y = [real(get(rx_sig,'FIELDY')-get(rx_sig,'FIELDY_TX')).';... 
             imag(get(rx_sig,'FIELDY')-get(rx_sig,'FIELDY_TX')).'];         
        f = [X;Y];
        
        return
        
end


%         function Essfm_backpropagation(dsp,Pavg,sig,Nspan)
%             r     = 10000/dsp.dz;           
%             
%             omega = 2*pi*sig.SYMBOLRATE*sig.FN'*1e9;             % [rad/s]
%             betaz = gpuArray(complex((0.5 * omega.^2*dsp.ch.b2...
%                                           + omega.^3*dsp.ch.b3/6)*dsp.dz));
%                                        
%             dx    = (dsp.ch.Lf>=dsp.dz)*dsp.dz     +... 
%                     (1-(dsp.ch.Lf>=dsp.dz))*dsp.ch.Lf;
%             k     = (dsp.nstep<1)*dsp.dz/dsp.ch.Lf +...
%                     (1-(dsp.nstep<1))*1;
%             
%             Leff  = dsp.dz;
%             if abs(dsp.ch.alphalin*dsp.dz) > 1e-6
%                 Leff = -k*(1-exp(dsp.ch.alphalin*dx))/dsp.ch.alphalin;
%             end
%             
%             z     = dsp.dz *(0:ceil(dsp.nstep)-1);
%             xi    = gpuArray(-dsp.ch.gamma*Leff*exp(dsp.ch.alphalin*z)*Pavg);
%             
% 
%             fiberpartion = (dsp.nstep<=1)*1+ (1-(dsp.nstep<=1))*Nspan;
%             
%             set(sig,'FIELDX', gpuArray(complex(get(sig,'FIELDX'))));
%             set(sig,'FIELDY', gpuArray(complex(get(sig,'FIELDY'))));
%             
%             for i=1:fiberpartion
%                 essfm(dsp,sig,xi,Nspan,betaz,r);
%             end
%             
%             set(sig,'FIELDX', gather(get(sig,'FIELDX')));
%             set(sig,'FIELDY', gather(get(sig,'FIELDY')));
%         
%         end
