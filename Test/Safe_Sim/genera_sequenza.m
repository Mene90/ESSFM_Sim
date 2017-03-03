function a = genera_sequenza(modul,pattern,dbl)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
usrad2=1./sqrt(2.);
switch modul
    case 'bpsk',
        if (strcmp(pattern,'debruijn'))
            % Sequenza di de Bruijn di '1' e '-1' di lunghezza 2^dbl:
            if (dbl<16)
                as=(idinput(2^dbl-1,'prbs'))';
                q=strfind(as,ones(1,dbl-1));
                a=[as(1:q-1),1,as(q:end)];
            elseif (dbl==16)
                a=load('debruijn16');
            elseif (dbl==17)
                a=load('debruijn17');
            elseif (dbl==18)
                a=load('debruijn18');
            end
        elseif (strcmp(pattern,'random'))
            % Sequenza casuale di '1' e '-1'
            a=2*round(rand(1,dbl))-1;
        end
    case 'qpsk',
        if (strcmp(pattern,'debruijn'))
            % Sequenza di de Bruijn di '1' e '-1' di lunghezza 2^dbl2:
            dbl2=2*dbl;
            if (dbl2<16)
                as=(idinput(2^dbl2-1,'prbs'))';
                q=strfind(as,ones(1,dbl2-1));
                ar=[as(1:q-1),1,as(q:end)];
            elseif (dbl2==16)
                ar=load('debruijn16');
            elseif (dbl2==17)
                ar=load('debruijn17');
            elseif (dbl2==18)
                ar=load('debruijn18');
            end
            a=usrad2*(ar+1i*circshift(ar,[1,dbl]));
        elseif (strcmp(pattern,'random'))
            % Sequenza casuale di '1/sqrt(2)' e '-1/sqrt(2)' su entrambe le quadrature
            a=usrad2*((2*round(rand(1,dbl))-1)+1i*(2*round(rand(1,dbl))-1));
        end
    case 'pulse',
        a=zeros(1,dbl);
        a(ceil(dbl/2))=1;
    case 'WGN',
        % Sequenza di campioni Gaussiani circolari i.i.d. lunga dbl
        a=usrad2*(randn(1,dbl)+1i*randn(1,dbl));
        
end

return
end

