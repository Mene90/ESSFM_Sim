classdef ElectricSource
    %ELECTRICSOURCE Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        signal
        ptype
        duty
        roll
    end
    
    methods
        function obj = ElectricSource(signal, ptype,...
                duty, roll)
            if isa(signal, 'Signal')
                obj.signal  = signal;
            else
                error(['the first arg must be an instance of the',...
                    'class Signal']);
            end
                        
            if max(duty) > 1
                error('duty must be <= 1');
            end
            
            if strcmp(ptype,'cosroll') && (max(roll) > 1)
                error('roll-off must be <= 1');
            end
            

            obj.ptype   = ptype;
            obj.duty    = duty;
            obj.roll    = roll;
            
        end
        
        function [elec] = pat2electricsource(obj, pat, format)
            
%             if max(pat) > 1
%                 error('%s requires a binary pattern',format);
%             end
        
            switch format
                case {'ook'}
                    elec = obj.elecsrc(pat);
                case {'bpsk','qpsk'}
                    elec = obj.elecsrc(2*pat-1);
                otherwise
                    error('wrong modulation format');
            end
                     
        end
    end
    
    methods 
    
        function elec = elecsrc(obj, pat)
            NT      = obj.signal.NT;
            NSYMB   = obj.signal.NSYMB;
            Nfft    = NSYMB * NT;                       % number of points
            elec    = zeros(Nfft,1);
            elpulse = pulseshape(obj,NT);  % single pulse
            
            nstart  = Nfft-NT+1;    % The first elpulse starts ate the end 
                                    % of the sequence (cyclic periodicity)
                                    % CAMBIA SE NON HA SENSO
            nend    = Nfft;
            
            elec(nstart:nend)   = pat(1)*elpulse(1:NT);
            elec(1:NT)          = pat(1)*elpulse(NT+1:NT*2);
            
            for kbit=2:NSYMB          
                nstart = (kbit-2)*NT+1;
                nend = kbit*NT;
                elec(nstart:nend) = elec(nstart:nend)+pat(kbit)*elpulse; 
            end
            
        end
        
        function y=pulseshape(obj,Nt)
            
            %PULSESHAPE Creates the fundamental pulse
            % Y=PULSESHAPE(NT) returns in Y a vector
            % [2*NT,1] containing the fundamental pulse whose type is
            % defined in ptype. NT is the number of points per bit.
            % ROLL is the pulse roll-off,
            % The length of Y is 2*NT because the 
            % roll-off spreads the pulse outside the bit time.
            
            ptype = obj.ptype;
            roll  = obj.roll;
            duty  = obj.duty;
            
            elpulse = zeros(Nt*2,1);     % elementary pulse (over two 
                                         % bit times because the roll-off 
                                         % spreads the pulse).
            switch ptype
                case 'cosroll'
                    
                    nl = round(0.5*(1-roll)*duty*Nt); % start index of cos
                                                      % roll-off
                                                         
                    nr = duty*Nt-nl-1;                % end index of cos
                                                      % roll-off
                    
                    nmark = 1:nl;      % indexis where the pulse is 1
                                                                     
                    ncos  = nl:nr;     % transition region of cos roll-off
                    
                    elpulse(Nt+nmark) = 1;
                    hperiod = duty*Nt-2*nl;
                    
                    if hperiod ~= 0
                        elpulse(ncos+Nt+1) = 0.5*(1+cos(pi/(hperiod)*...
                            (ncos-nl+0.5)));
                    end
                    
                    elpulse(1:Nt) = flipud(elpulse(Nt+1:Nt*2)); % first half of the pulse
                    
                case 'rootcosroll'
                    
                    nl = round(0.5*(1-roll)*duty*Nt); % start index of cos
                                                      % roll-off
                                                         
                    nr = duty*Nt-nl-1;                % end index of cos
                                                      % roll-off
                    
                    nmark = 1:nl;      % indexis where the pulse is 1
                                                                     
                    ncos  = nl:nr;     % transition region of cos roll-off
                    
                    elpulse(Nt+nmark) = 1;
                    hperiod = duty*Nt-2*nl;
                    
                    if hperiod ~= 0
                        elpulse(ncos+Nt+1) = sqrt(0.5*(1+cos(pi/(hperiod)*...
                            (ncos-nl+0.5))));
                    end
                    
                    elpulse(1:Nt) = flipud(elpulse(Nt+1:Nt*2)); % first half of the pulse
                
                case 'dirac'
                    
                    elpulse(Nt+1) = 1;
                
                otherwise
                    
                    error('error in pulseshape.m: the pulse ptype does not exist');
            
            end
            
            y=elpulse;
        end
        
    end
    
end

