function [ data ] = BER_essfm_vs_ssfm(Nstep,NC,dBm,size_length,n_prop_steps,etasp)
%BER_essfm_vs_ssfm Confronta la BER ottenuta utilizzando l'algoritmo
%split-spet Fourier classico (SSFM) o migliorato (ESSFM) per backpropagation su link multispan.
%
% Prima viene ottimizzata la configurazione dello ESSFM minimizzando MSE su un blocco di
% dati, poi viene fatta una simulazione su diversi blocchi per stimare BER
% Questa funzione viene chiamata da un lanciatore nel quale vengono settati
% opportuni dati di ingresso

% Alcune costanti utili:
clight = 2.99792458e8;	 % [km/s]
hplank = 6.62606957e-34; % [J*s]
lambda = 1.55e-6;        % [um]

%% Dati di ingresso:
%
% Nome del file di output:
%output_name='ESSFM_L30x120_QPSK_Rs32_filtr7_NS05_derot'; (non usato, i
%risultati vengono passati al lanciatore che provvede a salvarli)

% Parametri del link:
alf=0.2;
al=alf*0.230258509299405*1e-3;        % attenuation parameter (1/m)
b2=-21.67e-27; %parametro di GVD
gm=1.27e-3; %parametro di nonlinearità
LL=1.2e5;      %lunghezza di propagazione
Nspan=40;      %numero di span identiche
% etasp=2;        %coeff. di emissione spontanea degli amplificatori

% Caratteristiche del segnale trasmesso
modul='qpsk';   %'bpsk', 'qpsk', 'WGN'. Usato per ottimizzazione
modul_BER='qpsk';   %'bpsk', 'qpsk', 'WGN'. Usato per calcolo BER
pattern='debruijn'; %'debruijn', 'random', 'single'
pattern_BER='random'; %'debruijn', 'random', 'single'. Usato per calcolo BER
Rs=32e9;    %symbol rate [Hz]
Ps_dBm = dBm;  %Range di potenze medie in ingresso da considerare [dBm]
Ps=10.^(0.1*(Ps_dBm-30));    %potenza media (energia di ciascun impulso base/intervallo di segnalazione) [W]
Plen=length(Ps);            %numero di valori di potenza da testare
pls.shape='RRC';   %tipo di impulso
pls.ord=0.2;       %ordine dell'impulso supergaussiano
pls.bw=1;     %banda dell'impulso normalizzata alla symbol rate B/Rs
Nxs=2;          %Numero di campioni per tempo di simbolo
%  E' necessario scegliere Nxs in modo che la trasformata dell'impulso sia (approssimativamente)
% nulla al di fuori della banda di frequenze considerata [-Nxs*Rs/2,Nxs*Rs/2].
dbl=5;        %Lunghezza base per ottimizzazione: BPSK con de Bruijn diventa 2^dbl, QPSK con de Bruijn diventa 2^(2*dbl), altrimenti rimane dbl
dbl_BER=size_length;      %Lunghezza base per calcolo BER: BPSK con de Bruijn diventa 2^dbl, QPSK con de Bruijn diventa 2^(2*dbl), altrimenti rimane dbl
% Parametri degli algoritmi SSFM e ESSFM:
Nstep_exact=n_prop_steps;      %Numero passi SSFM per span per la propagazione (tanti per avere soluzione "esatta")
%Nstep=0.5;       %Numero di passi SSFM e ESSFM per span per la
%backpropagation. Può essere minore di 1 per indicare passi maggiori di una
%span (es.: 0.5 vuol dire un passo ogni 2 span). Qui non è attivo perché
%passato dal lanciatore.
Nst_tot=Nstep*Nspan;
if (round(Nst_tot)~=Nst_tot)
    display('ATTENZIONE! numero totale di passi richiesti non intero!');
end
tipo_ESSFM='filtr'; %'filt','filtr','quadform','quadformr'
%NC=4;               %numero di coefficienti (non attivo, passato dal
%lanciatore)
% solo ottimizzazione globale

%Opzioni di ottimizzazione:
options = optimset('Algorithm','trust-region-reflective','Display','off','Jacobian','off','DerivativeCheck','off','TolFun',1e-13,'TolX',1e-13);


%% Calcolo preliminare grandezze e vettori utili per la propagazione (che non cambiano all'interno del loop sulle potenze)

% Sequenze da trasmettere:
a=genera_sequenza(modul,pattern,dbl);                   % usata per ottimizzazione
a_tx=genera_sequenza(modul_BER,pattern_BER,dbl_BER);    % usata per calcolo BER
ar_tx=(real(a_tx)>=0)*2-1;                              % simboli della componente in fase
ai_tx=(imag(a_tx)>=0)*2-1;                              % simboli della componente in quadratura

% Segnale modulato, normalizzato a potenza unitaria:
u=linpulses(a,Nxs,pls);
utx=linpulses(a_tx,Nxs,pls);
Tc=1.0/(Rs*Nxs);          %intervallo di campionamento selezionato

%Filtro matchato:
% - per ottimizzazione:
N=length(u);
df=Nxs/N;
N1=floor(N/2);
N2=ceil(N/2)-1;
f=df*[0:N2,-N1:-1];       %vettore delle frequenze normalizzate dopo la FFT
Hf=filt(pls,f);
if (N2>N1)
    Hf(N1+1)=0.;
end
% - per calcolo BER:
N=length(utx);
df=Nxs/N;
N1=floor(N/2);
N2=ceil(N/2)-1;
f=df*[0:N2,-N1:-1];       %vettore delle frequenze normalizzate dopo la FFT
Hf_BER=filt(pls,f);
if (N2>N1)
    Hf_BER(N1+1)=0.;
end

%Varianza del rumore ASE da aggiungere e rapporto segnale rumore complessivo:
Gm1=(exp(al*LL)-1.d0);                          %Gain of the amplifier minus 1 (G-1)
N0=Nspan*Gm1*hplank*clight/lambda*etasp;        %Accumulate ASE Noise PSD
SNR=Ps/Rs/N0;                                   %Es/N0
% SNRdB=10*log10(SNR);


%% Ciclo sui valori di potenza in ingresso
% err_fin_st=zeros(1,Plen);
% err_fin_enh=zeros(1,Plen);
% err_fin_st_derot=zeros(1,Plen);
% err_fin_enh_derot=zeros(1,Plen);
BER_st=zeros(1,Plen);
BER_enh=zeros(1,Plen);
Copt=zeros(NC,Plen);
tic
for nn=1:Plen
    
    %Segnale all'uscita della fibra con propagazione "esatta" (SSFM a tanti passi)
    urx=u;
    sgn=sqrt(Nxs*0.5/SNR(nn)/Nspan);                        %deviazione standard campioni di rumore per ogni quadratura
    for is=1:Nspan
        wn=sgn*(randn(1,length(urx))+1i*randn(1,length(urx)));  %rumore AWGN
        urx=urx+wn;
        urx=ssfm(urx,Ps(nn),Tc,alf,b2,gm,LL,Nstep_exact);                   %split step
    end

    
    %% Funzione di errore da minimizzare (usa un puntatore ad una funzione anonima per passare i parametri extra):
    %Funzioni di errore globale:  (!!!solo filtr è effettivamente ottimizzata con derotazione e filtro matchato, le  altre 3 vanno sistemate!)
    switch tipo_ESSFM
        case 'filt',
            fmin=@(C)essfm_nspan_filt_opt(urx,ufin,Ps(nn)*exp(-al*LL),Tc,-alf,-b2,-gm,LL,Nspan,Nstep,C);    %filtraggio dei moduli quadri con coeff. complessi
        case 'filtr',
            %              fmin=@(C)essfm_nspan_filtr_opt(urx,ufin,Ps(nn)*exp(-al*LL),Tc,-alf,-b2,-gm,LL,Nspan,Nstep,C);    %filtraggio dei moduli quadri con coeff. reali
            fmin=@(C)essfm_nspan_filtr_opt(urx,a,Nxs,Ps(nn)*exp(-al*LL),Tc,-alf,-b2,-gm,LL,Nspan,Nstep,C,Hf);    %filtraggio dei moduli quadri con coeff. reali
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
    Copt(:,nn) = C;
end
tic
%% Simulazione per calcolo della BER
parfor nn=1:Plen    
        
    %Segnale all'uscita della fibra con propagazione "esatta" (SSFM a tanti passi)
    urx=utx;
    sgn=sqrt(Nxs*0.5/SNR(nn)/Nspan); 
    for is=1:Nspan
        wn=sgn*(randn(1,length(urx))+1i*randn(1,length(urx)));  %rumore AWGN
        urx=urx+wn;
        urx=ssfm(urx,Ps(nn),Tc,alf,b2,gm,LL,Nstep_exact);       %split step
    end
        
    %Calcolo segnale backpropagato con i due algoritmi:
    ufin_st=urx;
    ufin_enh=urx;
    if Nstep>=1,
        for is=1:Nspan,
            ufin_st = ssfm (ufin_st,Ps(nn)*exp(-al*LL),Tc,-alf,-b2,-gm,LL,Nstep);
            ufin_enh= essfm(ufin_enh,Ps(nn)*exp(-al*LL),Tc,-alf,-b2,-gm,LL,Nstep,tipo_ESSFM,Copt(:,nn));
        end
    else
        for is=1:round(Nspan*Nstep),
            ufin_st =ssfm(ufin_st,Ps(nn)*exp(-al*LL),Tc,0.,-b2,-gm*(1.0-exp(al*LL))/(-al*LL),LL/Nstep,1);
            ufin_enh=essfm(ufin_enh,Ps(nn)*exp(-al*LL),Tc,0.,-b2,-gm*(1.0-exp(al*LL))/(-al*LL),LL/Nstep,1,tipo_ESSFM,Copt(:,nn));
        end
    end
    
    
    % Filtro matchato:
    ufin_st=ifft(fft(ufin_st).*Hf_BER);
    ufin_enh=ifft(fft(ufin_enh).*Hf_BER);
    
    % Derotazione:
    rot_st=angle(mean(ufin_st(1:Nxs:end).*conj(a_tx)));
    rot_enh=angle(mean(ufin_enh(1:Nxs:end).*conj(a_tx)));
    ufin_st_derot=ufin_st*exp(-1i*rot_st);
    ufin_enh_derot=ufin_enh*exp(-1i*rot_enh);
    
    % misura BER:
    ar_st=(real(ufin_st_derot(1:Nxs:end))>=0)*2-1;         % Received symbol (real part)
    ar_enh=(real(ufin_enh_derot(1:Nxs:end))>=0)*2-1;       % Received symbol (real part)
    ai_st=(imag(ufin_st_derot(1:Nxs:end))>=0)*2-1;         % Received symbol (imaginary part)
    ai_enh=(imag(ufin_enh_derot(1:Nxs:end))>=0)*2-1;       % Received symbol (imaginary part)
    
    BER_st(nn)=0.5*(mean(ar_st~=ar_tx)+mean(ai_st~=ai_tx));
    BER_enh(nn)=0.5*(mean(ar_enh~=ar_tx)+mean(ai_enh~=ai_tx)); 
    
    display(['Power[dBm] = '   , num2str(dBm(nn))    ,char(9)...
             'BER con SSFM = ' , num2str(BER_st(nn)) ,char(9)...
             'BER con ESSFM = ', num2str(BER_enh(nn))            ]);
   
end
toc
data = [BER_st;BER_enh]';
% data=[Ps_dBm;SNRdB;BER_st;BER_enh]';
% data = Copt';
% save(output_name,'data','-ascii') %qui non salva, salva nel lanciatore


end

