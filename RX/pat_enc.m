function [ pat,patmat ] = pat_enc( pat )
%PAT_ENCODER Summary of this function goes here
%   Detailed explanation goes here


        stars_t = pat2stars(pat);
        stars_r = -ones(size(stars_t));

        for ii=2:length(stars_r)
            stars_r(ii) = -conj(stars_t(ii)).*stars_r(ii-1);
        end

        [pat,patmat] = stars2pat(stars_r);
        
        patmat=~patmat;
        pat=3-pat;

end

