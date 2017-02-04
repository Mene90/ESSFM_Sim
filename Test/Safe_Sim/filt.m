function H = filt(pls,f)
%FILT Restituisce la funzione di trasferimento di un dato filtro
%
%   *** ARGOMENTI DELLA FUNZIONE ***
%  pls:     caratteristiche dell'impulso (pls.shape,pls.ord,pls.bw)
%    f:     vettore delle frequenze
%
% banda e frequenze sono normalizzate allo stesso modo


if strcmp(pls.shape,'G')
    %impulso supergaussiano di ordine pls.ord e di banda (passabasso a 3dB) pls.bw
    H=exp(-0.5*log(2)*((f./(0.5*pls.bw)).^(2*pls.ord)));
elseif strcmp(pls.shape,'RC')
    fn=f/pls.bw;    %frequenza normalizzata alla banda
    ro=pls.ord;     %roll-off
    H=ones(1,length(fn));
    ind=find(abs(fn)>0.5*(1-ro))&(abs(fn)<0.5*(1+ro));
    H(ind)=0.5*(1+cos(pi/ro*(abs(fn(ind))-0.5*(1-ro))));
    H(abs(fn)>=0.5*(1+ro))=0;
elseif strcmp(pls.shape,'RRC')
    fn=f/pls.bw;    %frequenza normalizzata alla banda
    ro=pls.ord;     %roll-off
    H=ones(1,length(fn));
    ind=find((abs(fn)>0.5*(1-ro))&(abs(fn)<0.5*(1+ro)));
    H(ind)=sqrt(0.5*(1+cos(pi/ro*(abs(fn(ind))-0.5*(1-ro)))));
    H(abs(fn)>=0.5*(1+ro))=0;
end

