        
        HPLANCK = 6.62606896e-34;       % Planck's constant [J*s]
        CLIGHT = 299792458;             % speed of light [m/s]

        
        link.Nspan       = 10;
        link.LL          = 1e5;
        link.attenuation = 0.2;
        link.lambda      = 1550;
        link.sprop       = 10;
        link.nlindex     = 0;%2.5e-20;
        link.disp        = 17;
        
        sp.bprop    = 1;%link.sprop;
        
        amp.type    = 'Raman';
        amp.etasp   = 1;
        amp.Fn      = 10*log10(2*amp.etasp);
        
        al =link.attenuation*0.230258509299405*1e-3;
        Gm1=(exp(al*link.LL)-1.d0);
        N0=link.Nspan*Gm1*HPLANCK*CLIGHT/(link.lambda* 1e-9)*amp.etasp;
        
        pdbm                 = (-13:1:-4);%[-15,-10,-5,-3,-1,0,1,3,5,10];
        signal_prop.nt       = 6;
        signal_prop.nc       = 5;
        signal_prop.nsymb    = 2^19;
        signal_prop.symbrate = 50;    
        
        wdm.cch              = 3;
        
        snr0dB = zeros(length(pdbm),1);
        
        for i=1:length(pdbm)
            [signals{i},snr0dB(i),ch] = Test_mux(link,sp,signal_prop,amp,pdbm(i),wdm);
            sgs(i).snr0dB = snr0dB;
            sgs(i).sg = 1;
            sgs(i).Pu = pdbm(i);
            sgs(i).subc(1).tx = signals{i}.FIELDX_TX;
            sgs(i).subc(1).rx = signals{i}.FIELDX;
        end
        
        ch_properties       = ch.getProperties;
        ch_properties.Nspan = link.Nspan;
        
        savefile = strcat('G','_',int2str(link.LL/1000),'X',int2str(link.Nspan),'_WDM_',int2str(signal_prop.nc),'_',amp.type,'_nt_',int2str(signal_prop.nt));

        save(savefile,'sgs','ch_properties','amp','signal_prop','pdbm');