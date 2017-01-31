classdef Channel
      
    properties (SetAccess = private)
        signal      % signal to be transmited
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
        
        nplates 
        db0 
        theta
        epsilon
        dgdrms
        dgd
        db1
        lcorr       % waveplate length [m] (PMD step)
    end
    
    properties (Constant)
        CLIGHT      = 299792458;        % speed of light in vacuum [m/s]
        DEF_PLATES  = 100;
        sig0        = eye(2);       
        sig2        = [0 1;1 0];
        sig3i       = [0 1;-1 0];       % = i*sig3 = i*[0 -i;i 0]
    end
    
    methods
        
        function ch = Channel(LL,alphadB,lambda,aeff,n2,dispersion,...
                slope,nstep,pmd,nplates)
                        
            ch.nstep    = nstep;                   % step of SSFM    
            
            ch.lambda   = lambda;                  % wavelength [nm]
            ch.Lf       = LL;                      % fiber length [m]
            ch.disp     = dispersion;              % dispersion [ps/nm/km]
            ch.slope    = slope;                   % slope [ps/nm^2/km]
            ch.dz       = LL/nstep;                % step size
            ch.n2       = n2;                      % nonlinear index
            ch.aeff     = aeff;                    % effective area [um^2]
            
            ch.alphalin = (log(10)*1e-4)*alphadB;               % [m^-1]
            
            ch.gamma    = 2*pi*ch.n2 /...
                (lambda * ch.aeff) * 1e21;                        % [1/W/m]
            
%             b2 = -21.67e-27;
             ch.b2   = -ch.lambda^2 /...
                  2/pi/ch.CLIGHT * ch.disp* 1e-24;               % [s^2/m]
             
             ch.b3   = 0; 
%              b3   = (obj.lambda /(2*pi*obj.CLIGHT))^2 *...
%                  (2*obj.lambda*obj.disp + ...
%                  obj.lambda^2*obj.slope) * 1e-6;                  % [s^3/m]

%             b3   = (obj.lambda /(2*pi*obj.CLIGHT))^2 *...
%                 (2*obj.lambda*obj.disp + ...
%                 obj.lambda^2*obj.slope) * 1e-6;                  % [ns^3/m]
%             
%             omega       = 2*pi*sig.SYMBOLRATE*sig.FN';           % [rad/ns]
                       
            if(pmd)
                ch.nplates     = nplates;
                ch.db0         = rand(nplates,1)*2*pi - pi;
                ch.theta       = rand(nplates,1)*pi - 0.5*pi;      % azimuth: uniform R.V.;
                ch.epsilon     = 0.5*asin(rand(nplates,1)*2-1);    % uniform R.V. over the Poincare sphere
                ch.dgdrms      = sqrt((3*pi)/8)*dgd/sqrt(nplates);  
                ch.dgd         = dgd;
                ch.dgdrms      = ch.dgdrms/ch.signal.SYMBOLRATE; % convert in [ns]
                ch.db1         = ch.dgdrms*omega; 
                ch.lcorr       = ch.Lf/ch.nplates;
            end
        end
        
        function sig = propagates(ch,Pavg,sig)            
            ch_str  = struct('b2',ch.b2,'b3',ch.b3,'dz',ch.dz,'alphalin',ch.alphalin,'gamma',ch.gamma,'nstep',ch.nstep);
            sig_str = struct('FIELDX',sig.FIELDX','SYMBOLRATE',sig.SYMBOLRATE,'FN',sig.FN);
            sig_str = f_scalar_ssfm_mex(ch_str,Pavg,sig_str);
            set(sig,'FIELDX',sig_str.FIELDX');
        end
  
        function sig = scalar_ssfm(ch,Pavg,sig)
            
            omega = 2*pi*sig.SYMBOLRATE*sig.FN'*1e9;         % [rad/s]
            beta  = 0.5*omega.^2*ch.b2 + omega.^3*ch.b3/6;
            
            ux    = get(sig,'FIELDX');
            
            if abs(ch.alphalin*ch.dz) > 1e-6
                Leff     = (1-exp(-ch.alphalin*ch.dz))/ch.alphalin;
            else
                Leff     = ch.dz;
            end
            
            z    = ch.dz*(0:ch.nstep-1);
            
            xi   = ch.gamma*Leff*exp(-ch.alphalin*z)*Pavg;
            
            halfdz      = ch.dz/2;
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %                           HALF DZ GVD                      %
            %                                DZ SPM                      %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
             ux = scalar_lin_step(ch,beta*halfdz,ux); 
             ux = scalar_nl_step(ch,ux,xi(1)); 
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            for i=2:ch.nstep    
            
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %                           DZ GVD                       %
                %                           DZ SPM                       %
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                 ux = scalar_lin_step(ch,beta*ch.dz,ux); 
                 ux = scalar_nl_step(ch,ux,xi(i)); 
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %                      LAST HALF DZ GVD                      %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
             ux = scalar_lin_step(ch,beta*halfdz,ux);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            set(sig,'FIELDX',ux);
        end
        
        function sig = vectorial_ssfm(ch,Pavg,sig)
            
            omega = 2*pi*sig.SYMBOLRATE*sig.FN'*1e9;         % [rad/s]
            beta  = 0.5*omega.^2*ch.b2 + omega.^3*ch.b3/6;
            
            ux    = get(sig,'FIELDX');
            uy    = get(sig,'FIELDY');
            
            if abs(ch.alphalin*ch.dz) > 1e-6
                Leff     = (1-exp(-ch.alphalin*ch.dz))/ch.alphalin;
            else
                Leff     = ch.dz;
            end
            
            z    = ch.dz*(0:ch.nstep-1);
            
            xi   = ch.gamma*Leff*exp(-ch.alphalin*z)*Pavg;
            
            halfdz      = ch.dz/2;
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %                           HALF DZ GVD                      %
            %                                DZ SPM                      %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            [ux,uy] = vec_lin_step(ch,beta*halfdz,ux,uy); 
            [ux,uy] = vec_nl_step(ch,xi(1),ux,uy); 
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            for i=2:ch.nstep 
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %                           DZ GVD                       %
                %                           DZ SPM                       %
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                [ux,uy] = vec_lin_step(ch,beta*ch.dz,ux,uy); 
                [ux,uy] = vec_nl_step(ch,xi(i),ux,uy); 
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %                      LAST HALF DZ GVD                      %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            [ux,uy] = vec_lin_step(ch,beta*halfdz,ux,uy);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            set(sig,'FIELDX',ux);
            set(sig,'FIELDY',uy);
            
        end
        
        function [ux,uy] = vectorial_ssfm_pmd(ux,uy)
                        
          
            gam =obj.gamma*8/9;                        % 8/9 Manakov factor
            dzb = waveplatexstep(obj.dz,obj.lcorr);
            ntot = 0;
            
            halfntrunk  = ceil(0.5*obj.dz/obj.lcorr);  % trunks in half dz
            ntrunk      = ceil(obj.dz/obj.lcorr);      % trunks in dz
            
            halfdz      = 0.5*obj.dz;
            halfalpha   = 0.5*obj.alphalin;
            zprop       = 0.5*halfdz; 
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %                       HALF DZ GVD + PMD                    %
            %                            DZ SPM                          %
            %                       HALF DZ LOSS                         %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            [ux,uy] = ssfm_ln_step_pmd(dzb(1:halfntrunk),ntot,ux,uy);
            [ux,uy] = ssfm_nl_step_pmd(gam,obj.Leff,ux,uy);
            ux = ux*exp(-halfalpha*halfdz);    
            uy = uy*exp(-halfalpha*halfdz);    
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            ntot  = ntot + halfntrunk;
            zprop = zprop+obj.dz;
            step = 2;
            
            while zprop < obj.Lf         
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %                        DZ GVD + PMD                    %
                %                        DZ SPM                          %
                %                        DZ LOSS                         %
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                [ux,uy] = ssfm_ln_step_pmd(dzb(step,:),ntot,ux,uy);
                [ux,uy] = ssfm_nl_step_pmd(gam,obj.Leff,ux,uy);
                ux = ux*exp(-halfalpha*obj.dz);    
                uy = uy*exp(-halfalpha*obj.dz);
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                %%%%%%%%%%%%%%%%%%%% UPDATE PARAMETERS %%%%%%%%%%%%%%%%%%%%
                ntot  = ntot + ntrunk - logical(dbz(step,1));      
                zprop = zprop+obj.dz;
                step  = step+1;
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %                    LAST HALF DZ GVD + PMD                  %
            %                    LAST       DZ SPM                       %
            %                    LAST HALF DZ LOSS                       %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            [ux,uy] = ssfm_ln_step_pmd(dzb(step:halfntrunk),ntot,ux,uy);            
            [ux,uy] = ssfm_nl_step_pmd(gam,obj.Leff,ux,uy);             
            ux = ux*exp(-halfalpha*last_step);    
            uy = uy*exp(-halfalpha*last_step);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
        end
        
    end
    
    methods (Access = private)
        
        function ux         = scalar_lin_step(ch,betaxdz,ux)
            
            Hf = exp(-1i*betaxdz);
            ux = ifft( fft(ux) .* Hf);
            
        end
        
        function ux         = scalar_nl_step(ch,ux,xi)
                        
            pow = real(ux).^2 + imag(ux).^2;
            ux  = ux.*exp(-1i*xi.*pow);
            
        end
        
        function [ux,uy]    = vec_lin_step(ch,betaxdz,ux,uy)
            
            Hf = exp(-1i*betaxdz);
            
            ux = ifft( fft(ux) .* Hf);
            uy = ifft( fft(uy) .* Hf);
            
        end
        
        function [ux,uy]    = vec_nl_step(ch,xi,ux,uy)
            
            pow = real(ux).^2 + imag(ux).^2 + real(uy).^2 + imag(uy).^2;
            ux = ux.*exp(-1i*xi.*pow);
            uy = uy.*exp(-1i*xi.*pow);
            
        end
        
        function [ux,uy]    = ssfm_ln_step_pmd(dzb,ntot,ux,uy)
            
            s0   = obj.sig0;
            s2   = obj.sig2;
            s3i  = obj.sig3i;
            th   = obj.theta;
            eps  = obj.epsilon;
            
            nmem = ~(dzb(1) == obj.lcorr);    % if the previous trunk is 
                                              % not finished then nmen will 
                                              % be equal to 1
            
            ux = fft(ux);
            uy = fft(uy);
            
            ntrunk = length(dbz);
            
            for k=1:ntrunk
                n = ntot+k-nmem;                 
                
                matRth      = cos(th(n))*s0 - sin(th(n))*s3i;    
                matRepsilon = complex(cos(eps(n))*s0, sin(eps(n))*s2);     
                
                matR = matRth*matRepsilon;       
                
                %   sig.FIELDX(k) is
                %   the electric field for the kth frequency, we have that
                %   matR*D*matR'*A is the linear PMD step, where D is the 
                %   diagonal matrix where the DGD operates.
             
                
                uux = conj(matR(1,1))*ux + conj(matR(2,1))*uy;
                uuy = conj(matR(1,2))*ux + conj(matR(2,2))*uy;
                
                % 2) apply birefringence, DGD and GVD: 
                %    all in a diagonal matrix
                
                combeta     = obj.beta*dzb(k);   
                deltabeta   = 0.5*(obj.db1+obj.db0(n))...
                    *dzb(k)/obj.lcorr;           
                
                uux = exp(-1i*(combeta+deltabeta)).*uux;
                uuy = exp(-1i*(combeta-deltabeta)).*uuy;
                
                ux = matR(1,1)*uux + matR(1,2)*uuy;
                uy = matR(2,1)*uux + matR(2,2)*uuy;
            
            end
            
            ux = ifft(ux);
            uy = ifft(uy);
        end
        
        function [ux,uy]    = ssfm_nl_step_pmd(gam,leff,ux,uy)
            
            [ux,uy] = ssfm_nl_step(gam,leff,ux,uy);
            
        end
                
        function dzb        = waveplatexstep(dz,lcorr)
            
            nstep         = obj.Lf/obj.dz;
            dz_ntrunk     = ceil(dz/lcorr);
            halfdz_ntrunk = ceil(dz/2/lcorr);
            dzb           = ones([nstep+1,dz_ntrunk]);
            zprop = 0.5*dz;
            
            dzlast = 0.5*dz-lcorr*(halfdz_ntrunk-1);
            dzb(1,1:halfdz_ntrunk) = [lcorr*ones(1,halfdz_ntrunk-1),dzlast];
            dzmiss = lcorr - dzlast;
            
            zprop = zprop + dz;
            nthstep  = 2;
            
            while zprop < obj.Lf
                
                lastrunk   = lcorr*logical((lcorr - dzlast)/lcorr);
                dzlast     = dz - dzmiss - lcorr*(dz_trunk-2) - lastrunk;
                
                if dzmiss == 0
                    dzb(nthstep,:) = [lcorr*ones(1,dz_ntrunk-1),dzlast];
                else
                    dzb(nthstep,:) = [dzmiss,lcorr*ones(1,dz_ntrunk-2),dzlast];
                end
                
                dzmiss    = lcorr - dzlast;
                nthstep   = nthstep + 1;
                zprop     = zprop+dz;
            
            end
            
            dzlast = 0.5*dz-dzmiss-lcorr*(halfdz_ntrunk-1);
            
            if dzmiss == 0
                dzb(nthstep,1:halfdz_ntrunk) = ...
                    [lcorr*ones(1,halfdz_ntrunk-1),dzlast];
            else
                dzb(nthstep,1:halfdz_ntrunk) = ...
                    [dzmiss,lcorr*ones(1,halfdz_ntrunk-2),dzlast];
            end
            
        end
        
    end
    
end