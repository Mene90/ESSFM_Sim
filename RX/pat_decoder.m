function [ pat, patmat ] = pat_decoder( pat )
%PAT_DECODER Summary of this function goes here
%   Detailed explanation goes here

       
        stars_t = pat2stars(pat);
        
        stars_r=conj(stars_t).*circshift(stars_t,1);
        
        [pat, patmat] = stars2pat(stars_r);
        
        patmat=~patmat;     % both pattern have to be inverted
        pat=3-pat;
        
end

