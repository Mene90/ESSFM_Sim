function [ pat_rx ] = samp2pat( phases )
%SAMP2PAT Summary of this function goes here
%   Detailed explanation goes here


        second_bit = phases > 0;
        first_bit = abs(phases) <= pi/2;
        if size(phases,2) == 1
            pat_rx = [first_bit second_bit];
        else
            pat_rx = [first_bit(:,1) second_bit(:,1)  first_bit(:,2)  second_bit(:,2) ];
        end


end

