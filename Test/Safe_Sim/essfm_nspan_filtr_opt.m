function f = essfm_nspan_filtr_opt(x,xout,Nxs,Px,Tc,alf,bt2,gm,L,Nspan,Ns,C,Hf)
%
% Funzione errore da minimizzare. Restituisce l'errore commesso per ogni
% campione. Viene chiamata dalla subroutine lsqnonlin, che minimizza la somma dei quadrati degli errori.
%
% L'errore è calcolato backpropagando il segnale su link multispan,
% filtrandolo col filtro matchato, correggendo un eventuale errore di fase
% costante residuo (dovuto a rumore di fase nonlineare medio) ed infine
% confrontando il risultato con i simboli trasmessi

a=alf*0.230258509299405*1e-3;        % attenuation parameter (1/m)
y=x;
if Ns>=1
    for is=1:Nspan,
        y=essfm(y,Px,Tc,alf,bt2,gm,L,Ns,'filtr',C);
    end
else
    for is=1:round(Nspan*Ns)
        y=essfm(y,Px,Tc,0.,bt2,gm*(1.0-exp(-a*L))/(a*L),L/Ns,1,'filtr',C);
    end
end

% Filtro matchato e campionamento
y=ifft(fft(y).*Hf);
ys=y(1:Nxs:end);

% Derotazione
rot=angle(mean(ys.*conj(xout)));
ys=ys*exp(-1i*rot);

% Funzione errore
f=[real(ys-xout);imag(ys-xout)];
    
return
end

