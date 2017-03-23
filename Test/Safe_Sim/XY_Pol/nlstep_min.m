function [y1,y2] = nlstep_min(C,x,y)
% funzione che propaga il segnale in un passo nonlineare di fibra per
% l'algoritmo ESSFM. Applica una rotazione di fase pari ad una versione 
% filtrata della potenza istantanea in ingresso. Il filtro è simmetrico.

M=length(x);
N=length(C);

% Alloca spazio per il segnale
% y1=zeros(1,M);
% y2=zeros(1,M);
% Calcola moduli quadri (con prolungamento periodico)
xx=abs(x).^2 + abs(y).^2;
xx=[xx(M-N+2:M);xx;xx(1:N-1)];   %periodicamente
%xx=[0;xx;0];                     %con zeri

% Duplica il vettore dei coefficienti per avere risposta simmetrica
if length(C)>1
CC=[flipud(C);C(2:end)];
else
CC = flipud(C);
end
% Filtra il vettore dei moduli quadri con i coefficienti del filtro dato:
theta=conv(xx,CC,'valid');

% Calcola il segnale in uscita:
y1=transpose(x.*exp(1i*8/9*theta));
y2=transpose(y.*exp(1i*8/9*theta));

return
end
