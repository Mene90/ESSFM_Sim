%nlindex = [1.1e-20,2.5e-20,5e-20];

HPLANCK = 6.62606896e-34;       % Planck's constant [J*s]
CLIGHT = 299792458;             % speed of light [m/s]

n = input('Enter a number: ');

switch n
    case 1
        
%         distribution = input('Enter a distribution: ');
        Nspan            = [80,90,100,110,120];  
        link.Nspan       = Nspan(1);
        link.LL          = 0.8e5;
        link.attenuation = 0.2;
        link.lambda      = 1550;
        link.sprop       = 5;
        
        sp.bprop    = 1/link.Nspan;
        
        amp.etasp   = 2.4;
        amp.Fn      = 10*log10(2*amp.etasp);
        
        al =link.attenuation*0.230258509299405*1e-3;
        Gm1=(exp(al*link.LL)-1.d0);
        N0=link.Nspan*Gm1*HPLANCK*CLIGHT/(link.lambda* 1e-9)*amp.etasp;

        signal_prop.nsymb    = 2^20;
        signal_prop.symbrate = 50;
                
%         SNR_dB   = (-10:10:60);
%         pdbm     = 30+10*log10(signal_prop.symbrate*10^9*N0*10.^(SNR_dB*0.1));

        pdbm    =   [-5,-3,0,5,10,15,25,35,45];
        P       =   10.^((pdbm-30)*0.1);
        SNR     =   P/signal_prop.symbrate/10^9/N0;
        C       =   log2(1+SNR);
        SNR_dB  =   10*log10(SNR);
        distribution = 'HG';
        for k = 1:2
            if k == 2
                distribution = 'G';
            end
            for j = 1:length(Nspan)
                link.Nspan = Nspan(j);
                sp.bprop    = 1/link.Nspan;
                for i=1:length(pdbm)
                    [signals{i},SNRdB{i},ch] = TestZeroDispHighPower(link,sp,signal_prop,amp,pdbm(i),distribution);
                end
                
                ch_properties       = ch.getProperties;
                ch_properties.Nspan = link.Nspan;
                
                savefile        = strcat(distribution,'_',int2str(link.LL/1000),'X',int2str(link.Nspan),'_NoDisp_DBP_SS');
                
                save(savefile,'signals','SNRdB','ch_properties','amp','signal_prop','pdbm');
            end
        end
        
    case 2
        
        link.Nspan       = 30;
        link.LL          = 1.2e5;
        link.attenuation = 0.2;
        link.lambda      = 1550;
        link.sprop       = 30;
        link.nlindex     = 2.5e-20;
        link.disp        = 16;
        
        sp.bprop    = link.sprop;
        
        amp.etasp   = 1.76;
        amp.Fn      = 10*log10(2*amp.etasp);
        
        al =link.attenuation*0.230258509299405*1e-3;
        Gm1=(exp(al*link.LL)-1.d0);
        N0=link.Nspan*Gm1*HPLANCK*CLIGHT/(link.lambda* 1e-9)*amp.etasp;
        
        pdbm                 = (-10:4:8);
        signal_prop.nsymb    = 2^18;
        signal_prop.symbrate = 14;     
        
        for i=1:length(pdbm)
            [signals{i},SNRdB{i},ch] = TestConfronto(link,sp,signal_prop,amp,pdbm(i),'true');
        end
        
        ch_properties       = ch.getProperties;
        ch_properties.Nspan = link.Nspan;
    
        savefile        = strcat('G','_',int2str(link.LL/1000),'X',int2str(link.Nspan),'_DBP_Esatta');

        save(savefile,'signals','SNRdB','ch_properties','amp','signal_prop','pdbm');
        
    case 3

        link.Nspan       = 10;
        link.LL          = 1e5;
        link.attenuation = 0.2;
        link.lambda      = 1550;
        link.sprop       = 10;
        
        sp.bprop    = 5;
        
        amp.etasp   = 2;
        amp.Fn      = 10*log10(2*amp.etasp);
        
        al =link.attenuation*0.230258509299405*1e-3;
        Gm1=(exp(al*link.LL)-1.d0);
        N0=link.Nspan*Gm1*HPLANCK*CLIGHT/(link.lambda* 1e-9)*amp.etasp;

        signal_prop.nsymb    = 2^18;
        signal_prop.symbrate = 50;
                
%         SNR_dB   = (-10:10:60);
%         pdbm     = 30+10*log10(signal_prop.symbrate*10^9*N0*10.^(SNR_dB*0.1));

        pdbm    =   (-10:2:10);
        P       =   10.^((pdbm-30)*0.1);
        SNR     =   P/signal_prop.symbrate/10^9/N0;
        C       =   log2(1+SNR);
        SNR_dB  =   10*log10(SNR);
        
        for i=1:length(pdbm)
             [signals{i},SNRdB{i},ch] = TestZeroDisp(link,sp,signal_prop,amp,pdbm(i));
        end
        
        ch_properties       = ch.getProperties;
        ch_properties.Nspan = link.Nspan;
    
        savefile        = strcat('ZeroDisp_',int2str(link.LL/1000),'X',int2str(link.Nspan),'_NoDBP');

        save(savefile,'signals','SNRdB','ch_properties','amp','signal_prop','pdbm');
        
end