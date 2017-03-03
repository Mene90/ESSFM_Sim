function y = nlstep_min(C,x)
% funzione che propaga il segnale in un passo nonlineare di fibra per
% l'algoritmo ESSFM. Applica una rotazione di fase pari ad una versione 
% filtrata della potenza istantanea in ingresso. Il filtro è simmetrico.

M=length(x);
N=length(C);

% Alloca spazio per il segnale
y=zeros(M,1);

% Calcola moduli quadri (con prolungamento periodico)
xx=abs(x).^2;
xx=[xx(M-N+2:M);xx;xx(1:N-1)];   %periodicamente
%xx=[0;xx;0];                     %con zeri

% Duplica il vettore dei coefficienti per avere risposta simmetrica
CC=[flipud(C);C(2:end)];

% Filtra il vettore dei moduli quadri con i coefficienti del filtro dato:
theta=conv(xx,CC,'valid');

% Calcola il segnale in uscita:
y=x.*exp(1i*theta);



return
end
