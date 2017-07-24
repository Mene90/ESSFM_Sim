classdef Channel
      
        properties
        nstep       % number of step per spam
        disp        % dispersion [ps/nm/km]
        slope       % slope [ps/nm^2/km]
        lambda      % wavelength [nm]
        n2          % nonlinear index
        alphalin    % fiber attenuation [m^-1]
        gamma       % nonlinear coefficient [1/mW/m]
        b2          % [s^2/m]
        b3          % [s^3/m]
        aeff        % effective area [um^2]
        Lf          % fiber length [m]
        dz          % step size [m]
        alphadB
        end
        
        properties (Constant)
            CLIGHT      = 299792458;        % speed of light in vacuum [m/s]
        end
        
    methods
        
        function ch = Channel(LL,alphadB,lambda,aeff,n2,dispersion,...
                slope,nstep)
                        
            ch.nstep    = nstep;                   % step of SSFM           
            ch.lambda   = lambda;                  % wavelength [nm]
            ch.Lf       = LL;                      % fiber span length [m]
            ch.disp     = dispersion;              % D dispersion [ps/nm/km]
            ch.slope    = slope;                   % S slope [ps/nm^2/km]
            ch.dz       = LL/nstep;                % size of single step
            ch.n2       = n2;                      % nonlinear index
            ch.aeff     = aeff;                    % effective area [um^2]
            ch.alphadB  = alphadB;                 % loss coefficient [dB] 
            
            ch.alphalin = (log(10)*1e-4)*alphadB;                 % [m^-1]
            
            ch.gamma    = 2*pi*ch.n2 /(lambda * ch.aeff) * 1e21;  % [1/W/m]
            
            ch.b2       = -ch.lambda^2 /...
                          2/pi/ch.CLIGHT * ch.disp* 1e-24;        % [s^2/m]
             
            ch.b3       = 0; 
            
        end
        
        function prop = getProperties(ch)
            prop.LL    = ch.Lf;
            prop.gm    = ch.gamma;
            prop.b2    = ch.b2;
            prop.alpha = ch.alphadB;
            prop.disp  = ch.disp;            
        end
        
        function sing_span_propagation(ch,sig,gpu)
            
            ssfm(ch,get(sig,'POWER'),sig,gpu);
            
        end
        
        function gpu_propagation(ch,Nspan,ampli,sig)
            
            set(sig,'FIELDX', gpuArray(complex(get(sig,'FIELDX'))));
            set(sig,'FIELDY', gpuArray(complex(get(sig,'FIELDY'))));
            
            for i = 1:Nspan
                AddNoise(ampli,sig);
                sing_span_propagation(ch,sig,'true');
            end
            
            set(sig,'FIELDX', gather(get(sig,'FIELDX')));
            set(sig,'FIELDY', gather(get(sig,'FIELDY')));
            
        end
        
        function propagation(ch,Nspan,ampli,sig)
            for i = 1:Nspan
                AddNoise(ampli,sig);
                sing_span_propagation(ch,sig,'false');
            end
        end
        
    end
    
    methods 
        
        function ssfm(ch,Pavg,sig,gpu)
            
            omega  = 2*pi*sig.SYMBOLRATE*sig.FN'*1e9;         % [rad/s]   
            
            betaz  = complex((0.5 * omega.^2*ch.b2 + omega.^3*ch.b3/6)*ch.dz);                 
            betahz = complex((0.5 * omega.^2*ch.b2 + omega.^3*ch.b3/6)*ch.dz/2);
            
            geff   = ch.gamma*Pavg*ch.dz;
            if abs(ch.alphalin*ch.dz) > 1e-6
                geff = geff*(1-exp(-ch.alphalin*ch.dz))/ch.alphalin/ch.dz;
            end            
            z    = ch.dz*(0:ch.nstep-1);            
            xi   = geff*exp(-ch.alphalin*z);
           
            if(any(sig.FIELDY))
                xi = xi*8/9;
            end
            
            if gpu
               betaz  = gpuArray(betaz);
               betahz = gpuArray(betahz);
               xi     = gpuArray(xi);
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %                           HALF DZ GVD                      %
            %                                DZ SPM                      %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            lin_step(ch,betahz,sig); 
            nl_step(ch,xi(1),sig); 
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            for i=2:ch.nstep                 
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %                           DZ GVD                       %
                %                           DZ SPM                       %
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                lin_step(ch,betaz,sig); 
                nl_step(ch,xi(i),sig); 
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %                      LAST HALF DZ GVD                      %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            lin_step(ch,betahz,sig);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        end
        
        function lin_step(~,betaxdz,sig)           
            set(sig,'FIELDX',ifft( fft(sig.FIELDX).* exp(-1i*betaxdz)));
            set(sig,'FIELDY',ifft( fft(sig.FIELDY).* exp(-1i*betaxdz)));            
        end
        
        function nl_step(~,xi,sig)             
            pow = abs(sig.FIELDX).^2+abs(sig.FIELDY).^2;
            set(sig,'FIELDX', sig.FIELDX.*exp(-1i*xi.*pow));
            set(sig,'FIELDY', sig.FIELDY.*exp(-1i*xi.*pow));             
        end
    end
end