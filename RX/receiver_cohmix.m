classdef receiver_cohmix
    %RECEIVER_COHMIX Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        LO_Detuning = 0;
        LO_PhaseNoise;
        LO_Ecw = 1;
    end
    
    methods
        
        function obj = receiver_cohmix(LO_Ecw,LO_PhaseNoise,LO_Detuning)
            obj.LO_Ecw = LO_Ecw;
            obj.LO_PhaseNoise = LO_PhaseNoise;
            obj.LO_Detuning = LO_Detuning;
        end
        
        function Irx = receive(obj,ch,ofilt,efilt)
            
            Nfft    = length(ch.signal.FN);
            npoints = 1:Nfft;
            nind    = npoints';
            
            OBPF    = myfilter(ofilt.type,ch.signal.FN,0.5*ofilt.obw,ofilt.oord);
            
            [sigx,sigy] = obj.getfilteredSignal(OBPF,nind,ch.signal);
            
            %% local oscilator %%
            Elo = obj.LO_Ecw * exp(1j*(obj.LO_Detuning + obj.LO_PhaseNoise));
            
            % lowpass filter
            Hf = myfilter(efilt.type,ch.signal.FN,efilt.ebw,efilt.eord);
            
%              EmixOutX= 0.5*[ sigx  -     Elo         ...
%                              sigx  +     Elo          ...
%                              sigx  - 1j *Elo         ...
%                              sigx  + 1j *Elo ];
                         
                         
            EmixOutX = [sigx * 1j + Elo * 1j  ...
                    sigx     - Elo          ...
                    sigx * 1j - Elo         ...
                    -sigx     + Elo * 1j];
                         
            IricX = real( EmixOutX .* conj(EmixOutX) );
            IricX = [ IricX(:,1) - IricX(:,2) , IricX(:,3) - IricX(:,4) ];
            IricX = real(ifft(fft(IricX) .* Hf));

            Irx = IricX;
            
            if( ch.signal.FIELDY ~= 0)
                %
                %             EmixOutY = [sigy .* 1j + Elo .* 1j ...
                %                     sigy      - Elo         ...
                %                     sigy .* 1j - Elo        ...
                %                     -sigy      + Elo .* 1j    ];
                
                
                EmixOutY= 0.5*[ sigy  -     Elo         ...
                    sigy  +     Elo          ...
                    sigy  - 1j *Elo         ...
                    sigy  + 1j *Elo ];
                
                IricY = real( EmixOutY .* conj(EmixOutY) );
                IricY = [ IricY(:,1) - IricY(:,2) , IricY(:,3) - IricY(:,4) ];
                IricY = real(ifft(fft(IricY) .* Hf));
                Irx = [Irx IricY];
                
            end
            
        end
    end
    
    methods (Access = private)
        function [sigx,sigy] = getfilteredSignal(obj,OBPF,nind,sig)
            sigx = fft(sig.FIELDX);
            sigx = sigx(nind);
            sigx = sigx .* OBPF;
            sigx = ifft(sigx);
            if( sig.FIELDY ~= 0)
                sigy = fft(sig.FIELDY);
                sigy = sigy(nind);
                sigy = sigy .* OBPF;
                sigy = ifft(sigy);
            else
                sigy = 0;
            end
        end
    end
    
end


