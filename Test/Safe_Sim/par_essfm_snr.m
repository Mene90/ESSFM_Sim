function [ avg_snr ] = par_essfm_snr( Nstep,NC,dBm,size_length,n_prop_steps,etasp,gamma,R,Nspan )

% Alcune costanti utili:
clight = 2.99792458e8;	 % [km/s]
hplank = 6.62606957e-34; % [J*s]
lambda = 1.55e-6;        % [um]


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                         Link parameters                                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
alf=0.2;
al=alf*0.230258509299405*1e-3;        % attenuation           [1/m]
b2=-21.67e-27;                        % GVD parameter         [s^2/m]
gm=gamma;                             % nonlinear coefficient [1/W/m] �
LL=1e5;                               % sinlge span length    [m]
% Nspan=Nspan;                        % Number identical of spans
% etasp=2;                            % emission coeff.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                         Signal characteristics                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
modul='qpsk';              %'bpsk', 'qpsk', 'WGN'. Usato per ottimizzazione
modul_BER='qpsk';          %'bpsk', 'qpsk', 'WGN'. Usato per calcolo BER
pattern='debruijn';        %'debruijn', 'random', 'single'
pattern_BER='random';      %'debruijn', 'random', 'single'. Usato per calcolo BER
Rs=R;%32e9;                %symbol rate [Hz]
Ps_dBm = dBm;              %Range di potenze medie in ingresso da considerare [dBm]
Ps=10.^(0.1*(Ps_dBm-30));  %potenza media (energia di ciascun impulso base/intervallo di segnalazione) [W]
Plen=length(Ps);           %numero di valori di potenza da testare
pls.shape='RRC';           %tipo di impulso
pls.ord=0.1;               %ordine dell'impulso supergaussiano
pls.bw=1;                  %banda dell'impulso normalizzata alla symbol rate B/Rs
Nxs=2;                     %Numero di campioni per tempo di simbolo
%  E' necessario scegliere Nxs in modo che la trasformata dell'impulso sia (approssimativamente)
% nulla al di fuori della banda di frequenze considerata [-Nxs*Rs/2,Nxs*Rs/2].
dbl=5;                     %Lunghezza base per ottimizzazione: BPSK con de Bruijn diventa 2^dbl, QPSK con de Bruijn diventa 2^(2*dbl), altrimenti rimane dbl
dbl_BER=size_length;       %Lunghezza base per calcolo BER: BPSK con de Bruijn diventa 2^dbl, QPSK con de Bruijn diventa 2^(2*dbl), altrimenti rimane dbl
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                    ESSFM and SSFM parameters                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Nstep_exact=n_prop_steps;      %Numero passi SSFM per span per la propagazione (tanti per avere soluzione "esatta")
%Nstep=0.5;                    %Numero di passi SSFM e ESSFM per span per la
%backpropagation. Può essere minore di 1 per indicare passi maggiori di una
%span (es.: 0.5 vuol dire un passo ogni 2 span). Qui non è attivo perché
%passato dal lanciatore.

Nst_tot=Nstep*Nspan;
if (round(Nst_tot)~=Nst_tot)
    display('ATTENZIONE! numero totale di passi richiesti non intero!');
end

tipo_ESSFM='filtr';           %'filt','filtr','quadform','quadformr'
%NC=4;                        %numero di coefficienti (non attivo, passato dal lanciatore) solo ottimizzazione globale

%Opzioni di ottimizzazione:
options = optimset('Algorithm','trust-region-reflective','Display','off','Jacobian','off','DerivativeCheck','off','TolFun',1e-13,'TolX',1e-13);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Calcolo preliminare grandezze e vettori utili per la propagazione 

% Sequenze da trasmettere:
[patx(:,1), patmatx]    = Pattern.debruijn(1,4,2^10);
[paty(:,1), patmaty]    = Pattern.debruijn(8,4,2^10);
ax = transpose(1./sqrt(2.)*((2*patmatx(:,1)-1)+1i*(2.*patmatx(:,2)-1)));      % usata per ottimizzazione
ay = transpose(1./sqrt(2.)*((2*patmaty(:,1)-1)+1i*(2.*patmaty(:,2)-1)));      % usata per ottimizzazione

ax_tx  = genera_sequenza(modul_BER,pattern_BER,dbl_BER);    % usata per calcolo BER (pol x)
ay_tx  = genera_sequenza(modul_BER,pattern_BER,dbl_BER);    % usata per calcolo BER (pol y)
arx_tx = (real(ax_tx)>=0)*2-1;                              % simboli della componente in fase       (pol x)
aix_tx = (imag(ax_tx)>=0)*2-1;                              % simboli della componente in quadratura (pol y)
ary_tx = (real(ay_tx)>=0)*2-1;                              % simboli della componente in fase       (pol x)
aiy_tx = (imag(ay_tx)>=0)*2-1;                              % simboli della componente in quadratura (pol y)

% Segnale modulato, normalizzato a potenza unitaria:
ux  = linpulses(ax,Nxs,pls);
uy  = linpulses(ay,Nxs,pls);

utx = linpulses(ax_tx,Nxs,pls);
uty = linpulses(ay_tx,Nxs,pls);

Tc=1.0/(Rs*Nxs);                  % intervallo di campionamento selezionato

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                   Optimization Matched filter                          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N  = length(ux);
df = Nxs/N;
N1 = floor(N/2);
N2 = ceil(N/2)-1;
f  = df*[0:N2,-N1:-1];       
Hf = filt(pls,f);
if (N2>N1)
    Hf(N1+1)=0.;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                       BER Matched filter                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N      = length(utx);
df     = Nxs/N;
N1     = floor(N/2);
N2     = ceil(N/2)-1;
f      = df*[0:N2,-N1:-1];   
Hf_BER = gpuArray(filt(pls,f));
if (N2>N1)
    Hf_BER(N1+1)=0.;
end%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                      Amplifier Noise parameters                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Gm1 = (exp(al*LL)-1.d0);                       %Gain of the amplifier minus 1 (G-1)
N0  = Nspan*Gm1*hplank*clight/lambda*etasp;    %Accumulate ASE Noise PSD
SNR = Ps/Rs/N0;                                %Es/N0
% SNRdB=10*log10(SNR);

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                               TRAINING                                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% err_fin_st=zeros(1,Plen);
% err_fin_enh=zeros(1,Plen);
% err_fin_st_derot=zeros(1,Plen);
% err_fin_enh_derot=zeros(1,Plen);
BER_st=zeros(1,Plen);
BER_enh=zeros(1,Plen);
% Copt=zeros(NC,Plen);

urx=gpuArray(ux);
ury=gpuArray(uy);
sgn=sqrt(Nxs*0.5/SNR/Nspan);                        %deviazione standard campioni di rumore per ogni quadratura
for is = 1:Nspan
    wnx = sgn*(randn(1,length(urx))+1i*randn(1,length(urx)));  %rumore AWGN
    wny = sgn*(randn(1,length(ury))+1i*randn(1,length(ury)));
    urx = urx+wnx;
    ury = ury+wny;
    [urx,ury] = ssfm(urx,ury,Ps,Tc,alf,b2,gm,LL,Nstep_exact);                   %split step
end


%% Funzione di errore da minimizzare (usa un puntatore ad una funzione anonima per passare i parametri extra):
%Funzioni di errore globale:  (!!!solo filtr è effettivamente ottimizzata con derotazione e filtro matchato, le  altre 3 vanno sistemate!)
switch tipo_ESSFM
    case 'filt',
        fmin=@(C)essfm_nspan_filt_opt(urx,ufin,Ps(nn)*exp(-al*LL),Tc,-alf,-b2,-gm,LL,Nspan,Nstep,C);    %filtraggio dei moduli quadri con coeff. complessi
    case 'filtr',
        %              fmin=@(C)essfm_nspan_filtr_opt(urx,ufin,Ps(nn)*exp(-al*LL),Tc,-alf,-b2,-gm,LL,Nspan,Nstep,C);    %filtraggio dei moduli quadri con coeff. reali
        fmin=@(C)essfm_nspan_filtr_opt(urx,ury,ax,ay,Nxs,Ps*exp(-al*LL),Tc,-alf,-b2,-gm,LL,Nspan,Nstep,C,Hf);    %filtraggio dei moduli quadri con coeff. reali
    case 'quadform',
        fmin=@(C)essfm_nspan_opt(urx,ufin,Ps(nn)*exp(-al*LL),Tc,-alf,-b2,-gm,LL,Nspan,Nstep,C);    %fase nonlineare pari a forma quadratica dei campioni
    case 'quadformr',
        fmin=@(C)essfm_nspanr_opt(urx,ufin,Ps(nn)*exp(-al*LL),Tc,-alf,-b2,-gm,LL,Nspan,Nstep,C);    %fase nonlineare pari a forma quadratica dei campioni
end

%Inizializzazione dei coefficienti:
switch tipo_ESSFM
    case {'filt','filtr'},
        C0=zeros(NC,1);
        C0(1,1)=1;%-gm*Leff*Ps(nn);
    case {'quadform','quadformr'}
        NCr=floor(0.5*sqrt(NC));
        CC0=zeros(2*NCr+1);
        CC0(NCr+1,NCr+1)=1;%-gm*Leff*Ps(nn);
        C0=reshape(CC0,1,NC);
end
%Ottimizzazione dei coefficienti (nonlinear LSQ)
[C,err]=lsqnonlin(fmin,C0,[],[],options);
Copt(:,1) = gpuArray(C);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copt
%% Simulazione per calcolo della BER
% for nn=1:Plen    
        
    %Segnale all'uscita della fibra con propagazione "esatta" (SSFM a tanti passi)
    urx=gpuArray(utx);
    ury=gpuArray(uty);
    
    sgn=sqrt(Nxs*0.5/SNR/Nspan); 
    for is=1:Nspan
        wnx=sgn*(randn(1,length(urx))+1i*randn(1,length(urx)));  %rumore AWGN
        wny=sgn*(randn(1,length(ury))+1i*randn(1,length(ury)));  %rumore AWGN
        urx=urx+wnx;
        ury=ury+wny;
        [urx,ury]=ssfm(urx,ury,Ps,Tc,alf,b2,gm,LL,Nstep_exact);       %split step
    end
        
    %Calcolo segnale backpropagato con i due algoritmi:
    ufinx_st=urx;
    ufiny_st=ury;
    ufinx_enh=urx;
    ufiny_enh=ury;
    if Nstep>=1,
        for is=1:Nspan,
%             [ufinx_st,ufiny_st]   = ssfm (ufinx_st,ufiny_st,Ps(nn)*exp(-al*LL),Tc,-alf,-b2,-gm,LL,Nstep);
            [ufinx_enh,ufiny_enh] = essfm(ufinx_enh,ufiny_enh,Ps*exp(-al*LL),Tc,-alf,-b2,-gm,LL,Nstep,tipo_ESSFM,Copt(:,1));
        end
    else
        for is=1:round(Nspan*Nstep),
%             [ufinx_st,ufiny_st] =ssfm(ufinx_st,ufiny_st,Ps(nn)*exp(-al*LL),Tc,0.,-b2,-gm*(1.0-exp(al*LL))/(-al*LL),LL/Nstep,1);
            [ufinx_enh,ufiny_enh]=essfm(ufinx_enh,ufiny_enh,Ps*exp(-al*LL),Tc,0.,-b2,-gm*(1.0-exp(al*LL))/(-al*LL),LL/Nstep,1,tipo_ESSFM,Copt(:,1));
        end
    end
    
    
    % Filtro matchato:
    ufinx_st=ifft(fft(ufinx_st).*Hf_BER);
    ufiny_st=ifft(fft(ufiny_st).*Hf_BER);
    ufinx_enh=ifft(fft(ufinx_enh).*Hf_BER);
    ufiny_enh=ifft(fft(ufiny_enh).*Hf_BER);
    
    % Derotazione:
    rotx_st=angle(mean(ufinx_st(1:Nxs:end).*conj(ax_tx)));
    roty_st=angle(mean(ufiny_st(1:Nxs:end).*conj(ay_tx)));
    rotx_enh=angle(mean(ufinx_enh(1:Nxs:end).*conj(ax_tx)));
    roty_enh=angle(mean(ufiny_enh(1:Nxs:end).*conj(ay_tx)));
    ufinx_st_derot=ufinx_st*exp(-1i*rotx_st);
    ufiny_st_derot=ufiny_st*exp(-1i*roty_st);
    ufinx_enh_derot=ufinx_enh*exp(-1i*rotx_enh);
    ufiny_enh_derot=ufiny_enh*exp(-1i*roty_enh);
    
    snr_x = mySNR(gather(ufinx_enh_derot(1:Nxs:end)),gather(ax_tx));
    snr_y = mySNR(gather(ufiny_enh_derot(1:Nxs:end)),gather(ay_tx));

    avg_snr = 10*log10(0.5*(snr_x+snr_y));
    
%     % misura BER:
%     arx_st  =(real(ufinx_st_derot(1:Nxs:end))>=0)*2-1;         % Received symbol (real part)
%     ary_st  =(real(ufiny_st_derot(1:Nxs:end))>=0)*2-1;
%     arx_enh =(real(ufinx_enh_derot(1:Nxs:end))>=0)*2-1;         % Received symbol (real part)
%     ary_enh =(real(ufiny_enh_derot(1:Nxs:end))>=0)*2-1;         % Received symbol (real part)
%     aix_st  =(imag(ufinx_st_derot(1:Nxs:end))>=0)*2-1;         % Received symbol (imaginary part)
%     aiy_st  =(imag(ufiny_st_derot(1:Nxs:end))>=0)*2-1;
%     aix_enh =(imag(ufinx_enh_derot(1:Nxs:end))>=0)*2-1;         % Received symbol (imaginary part)
%     aiy_enh =(imag(ufiny_enh_derot(1:Nxs:end))>=0)*2-1;
%     
%     BER_st(nn) = 0.5*(0.5*(mean(arx_st~=arx_tx)+ mean(aix_st~=aix_tx))   + 0.5*(mean(ary_st~=ary_tx)  + mean(aiy_st~=aiy_tx)));
%     BER_enh(nn)= 0.5*(0.5*(mean(arx_enh~=arx_tx)+mean(aix_enh~=aix_tx)) + 0.5*(mean(ary_enh~=ary_tx) + mean(aiy_enh~=aiy_tx))); 
    
%     data(nn,:) = avgber;
%     display(['Power[dBm] = '   , num2str(dBm(nn))    ,char(9)...
%              'BER con SSFM = ' , num2str(BER_st(nn)) ,char(9)...
%              'BER con ESSFM = ', num2str(BER_enh(nn))            ]);
   
% end

%  data = [BER_st;BER_enh]';

end

