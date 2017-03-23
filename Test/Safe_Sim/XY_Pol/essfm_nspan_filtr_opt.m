function f = essfm_nspan_filtr_opt(xx,xy,xout,yout,Nxs,Px,Tc,alf,bt2,gm,L,Nspan,Ns,C,Hf)
%
% Funzione errore da minimizzare. Restituisce l'errore commesso per ogni
% campione. Viene chiamata dalla subroutine lsqnonlin, che minimizza la somma dei quadrati degli errori.
%
% L'errore è calcolato backpropagando il segnale su link multispan,
% filtrandolo col filtro matchato, correggendo un eventuale errore di fase
% costante residuo (dovuto a rumore di fase nonlineare medio) ed infine
% confrontando il risultato con i simboli trasmessi

a  = alf*0.230258509299405*1e-3;        % attenuation parameter (1/m)
y1 = xx;
y2 = xy;
if Ns>=1
    for is=1:Nspan,
        [y1, y2] = essfm(y1,y2,Px,Tc,alf,bt2,gm,L,Ns,'filtr',C);
    end
else
    for is=1:round(Nspan*Ns)
        [y1, y2] = essfm(y1,y2,Px,Tc,0.,bt2,gm*(1.0-exp(-a*L))/(a*L),L/Ns,1,'filtr',C);
    end
end

% Filtro matchato e campionamento
y1=ifft(fft(y1).*Hf);
y2=ifft(fft(y2).*Hf);

y1s=y1(1:Nxs:end);
y2s=y2(1:Nxs:end);

% Derotazione
rotx = angle(mean(y1s.*conj(xout)));
roty = angle(mean(y2s.*conj(yout)));

y1s = y1s*exp(-1i*rotx);
y2s = y2s*exp(-1i*roty);

% Funzione errore
f=[gather(real(y1s-xout));gather(imag(y1s-xout));gather(real(y2s-yout));gather(imag(y2s-yout))];
    
return
end

