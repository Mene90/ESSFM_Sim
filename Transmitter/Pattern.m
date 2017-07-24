classdef Pattern
    %PATTERN Summary of this class goes here
    %   Detailed explanation goes here
        
    methods (Static = true)
             
        function [pat,bmap] = debruijn(nseed,alphabet_size,nsymbols)
            q = log2(alphabet_size);
            
            if rem(q,1)
                error(['De Bruijn sequence can be implemented only ',...
                    ' for power of 2 alphabets']);
            end
            
            maxseed = nsymbols*(nsymbols-1)/4;              % largest seed
            
            if nseed(1) > maxseed
                error('nseed must be <= %d',maxseed); 
            end
            
            lg2 = log2(nsymbols);
            if rem(lg2,q)
                error(['The De Bruijn sequence does not exist! ',...
                    'log2(number of symbols)=%d is not a ',...
                    'multiple of log2(alphabet length)=%d.'],...
                    lg2,q);
            end
            [pat,bmap] = Pattern.debruijn_seq(log2(nsymbols)/q,nseed(1),q);
            bmap=bmap';
            
        end
                
        function [cmap,pat,bmap] = random(alphabet_size,nsymbols)        
            pat  = floor(alphabet_size.*rand(1,nsymbols));
            bmap = Pattern.mydec2bin(pat,log2(alphabet_size));
            tmp =  2*bmap-1;
            cmap = 1./sqrt(2.)*(tmp(:,1)+1i*tmp(:,2));
        end
        
        function [pat]      = gaussian(nsymbols)
            usrad2=1./sqrt(2.);
            pat = usrad2*(randn(nsymbols,1)+1i*randn(nsymbols,1));
        end
        
        function [pat]      = halfgaussian(nsymbols)
            hgpd  = makedist('HalfNormal');
            upd   = makedist('Uniform');
            pat   = random(hgpd,nsymbols,1).*exp(1i*random(upd,nsymbols,1)*2*pi);
        end
        
        function [pat]      = bin2complex(bmap)
            usrad2   = 1./sqrt(2.);
            pat      = usrad2*(bmap(:,1)+1i*bmap(:,2));
        end
                
    end
    
    methods (Static = true, Access = protected)
                
        function [y,tmat]=debruijn_seq(n,seed,q)
            
            %DEBRUIJN_SEQ De Bruijn sequence generator
            %   Y=DEBRUIJN_SEQ(N,SEED) generates a De Bruijn sequence of 2^N bits.
            %   A De Bruijn sequence contains all patterns of N bits exactly once.
            %
            %   SEED is the random seed (integer). For SEED > 2^(N-2) the resulting
            %   sequence is a pseudo-random shifted copy of a sequence with
            %   SEED <= 2^(N-2).
            %
            %   Reference:
            %
            %   [1] F. S. Annexstein, "Generating De Bruijn Sequences: An Efficient
            %   Implementation," IEEE Transaction on Computers, vol. 48, no. 2, Feb.
            %   1997.
            
            ns = q*n; % log2 length of the DeBruijn sequence
            N = 2^(ns-2);
            nseed = mod(seed,N);
            x = double(dec2bin(nseed,ns-2)-48); % convert to pattern
            y = double(Pattern.generate_debruijn(ns,x));
            
            
            if q > 1 % multilevel DeBruijn
                mvm = floor(ns/q*((1:q)-1));
                tmat = zeros(q,2^ns); % columns: binary DeBruijn
                for kk=1:q
                    tmat(kk,:) = circshift(y,[mvm(kk),0]);
                end
                y = 2.^(q-1:-1:0) * tmat ; % bin2dec conversion
            else
                tmat = y;         %if q=1 then bin representation of y is y
            end
            
            if seed > N % no more seeds -> apply a random delay shift
                nseed2 = ceil((seed-N+1)/N);
                nshift = mod(97*nseed2,2^ns-1)+1; % congruent random generator
                y = circshift(y,-nshift); % shift on the right
                tmat = circshift(tmat.',[-nshift 0]).';
            end
        end              
        
        function y=generate_debruijn(n,x)
            
            if n == 2
                y=[1,1,0,0] == 1;
            elseif n == 3
                if x(1) == 0
                    y = [1,0,1,1,1,0,0,0] == 1;
                else
                    y = [1,1,1,0,1,0,0,0] == 1;
                end
            else
                y= Pattern.next_debruijn(Pattern.generate_debruijn(n-1,x),2^(n-2)+(-1)^x(n-3),x(n-2));
            end
        end
        
        function y=next_debruijn(w,i,k)
            
            C = Pattern.xor_iscan(w);
            Cbar = not(C);
            part1 = C(1:i-k);
            part2 = Cbar(i+k:end);
            part3 = Cbar(1:i-1+k);
            part4 = C(i-k+1:end);
            y = [part1,part2,part3,part4];
        end
        
        function y=xor_iscan(x)
            
            % Being x=[x1,x2,...,xn] and y=[y1,y2,...,yn] it is:
            %       yk = x1 xor x2 xor ... xk
            
            z = cumsum(x);
            y = rem(z,2) == 1;
        end
        
        function y=mydec2bin(d,q)
            
            [f,e]=log2(max(d));
            y=rem(floor(d(:)*pow2(1-max(1,e):0)),2);
            cols = size(y,2);
            if q > cols
                y(:,q-cols+(1:cols)) = y(:,1:cols);
                y(:,1:(q-cols)) = zeros(size(y,1),q-cols);
            end
            
        end
        
    end
    
end

