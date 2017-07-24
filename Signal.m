classdef Signal < matlab.mixin.SetGet & matlab.mixin.Copyable
    
    properties
        FIELDX       % x-component of the electric field
        FIELDY       % y-component of the electric field
        FIELDX_TX    % x-component of the electric field (not propagated)
        FIELDY_TX    % y-component of the electric field (not propagated)
        FN           % normalized frequencies
        NSYMB        % number of symbols
        NT           % points x symbol
        NCH          % Number of Channels
        SYMBOLRATE   % The symbolrate will be initialized in ELECTRICSOURCE
        LAMBDA       % The range of wavelength used to propagate the field 
        POWER        % power of the signal
    end
    
    methods
        
        function sig = Signal(Nsymb,Nt,symbolrate,lambda,nch)
            stepf     = 1/Nsymb;
            sig.FN    = fftshift(-Nt/2:stepf:Nt/2-stepf); % normalized freq.
            sig.NSYMB = Nsymb;                            % number of symb.
            sig.NT    = Nt;                               % points x symbol
            sig.SYMBOLRATE = symbolrate;
            sig.FIELDX     = zeros(Nsymb*Nt,1);
            sig.FIELDY     = zeros(Nsymb*Nt,1);
            sig.FIELDX_TX  = zeros(Nsymb,1);
            sig.FIELDY_TX  = zeros(Nsymb,1);
            sig.LAMBDA     = lambda;
            sig.NCH        = nch;
        end
        
        function data = getproperties(sig)
            data.FIELDX    = get(sig,'FIELDX');
            data.FIELDX_TX = get(sig,'FIELDX_TX');
            if(any(sig.FIELDY))
                data.FIELDY    = get(sig,'FIELDY');
                data.FIELDY_TX = get(sig,'FIELDY_TX');
            end
            data.POWER     = get(sig,'POWER');
        end
                   
    end
    
end

