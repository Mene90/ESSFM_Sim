classdef Signal < matlab.mixin.SetGet & matlab.mixin.Copyable
    %SIGNAL Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        FIELDX       % x-component of the electric field
        FIELDY       % y-component of the electric field
        FIELDX_TX    % x-component of the electric field (not propagated)
        FIELDY_TX    % y-component of the electric field (not propagated)
        FN           % normalized frequencies
        NSYMB        % number of symbols
        NT           % points x symbol
        SYMBOLRATE   % The symbolrate will be initialized in ELECTRICSOURCE
        POWER        % power of the signal
    end
    
    methods
        
        function sig = Signal(Nsymb,Nt,symbolrate)
            stepf   = 1/Nsymb;
            sig.FN  = fftshift(-Nt/2:stepf:Nt/2-stepf); % normalized freq.
            sig.NSYMB = Nsymb;                          % number of symb.
            sig.NT    = Nt;                             % points x symbol
            sig.SYMBOLRATE = symbolrate;
        end
        
                   
    end
    
end

