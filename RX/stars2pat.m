function [pat, patmat] = stars2pat( stars )
%STARS2PAT Summary of this function goes here
%   Detailed explanation goes here

        pat(stars==1)=0;
        pat(stars==i)=1;
        pat(stars==-1)=3;
        pat(stars==-i)=2;
        
        patmat=zeros(length(pat),2);
        
        patmat(stars== 1,1)=0; patmat(stars== 1,2)=0;
        patmat(stars== i,1)=0; patmat(stars== i,2)=1;
        patmat(stars==-1,1)=1; patmat(stars==-1,2)=1;
        patmat(stars==-i,1)=1; patmat(stars==-i,2)=0;

end

