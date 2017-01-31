function [ stars ] = pat2stars( pat )
%PAT2STARS Summary of this function goes here
%   Detailed explanation goes here

        if max(max(pat))>1 || size(pat,2) ~= 2
            error('pattern must be a binary matrix with size [Nsymb,2]');
        end
        stars=zeros(length(pat),1);

        for ii=1:length(stars)
            if pat(ii,:)==[0 0] stars(ii)=1;end
            if pat(ii,:)==[0 1] stars(ii)=1i;end
            if pat(ii,:)==[1 1] stars(ii)=-1;end
            if pat(ii,:)==[1 0] stars(ii)=-1i;end
        end

end

