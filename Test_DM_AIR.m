        
        HPLANCK = 6.62606896e-34;       % Planck's constant [J*s]
        CLIGHT = 299792458;             % speed of light [m/s]

        
        link.Nspan       = 10;
        link.LL          = 1e5;
        link.attenuation = 0.2;
        link.lambda      = 1550;
        link.sprop       = 1000;
        link.nlindex     = 2.5e-20;
        link.disp        = 17;
        
        sp.bprop    = 100;
        
        amp.type    = 'Raman';
        amp.etasp   = 1;
        amp.Fn      = 10*log10(2*amp.etasp);
        
        al =link.attenuation*0.230258509299405*1e-3;
        Gm1=(exp(al*link.LL)-1.d0);
        N0=link.Nspan*Gm1*HPLANCK*CLIGHT/(link.lambda* 1e-9)*amp.etasp;
        
        nsc                  = [1,2,4];
        nt                   = [1,2,4];
        
        
        
        
        pdbm                 = (-13:1:-4);
        
        signal_prop.nt       = 8;
        signal_prop.nc       = 5;        
   
        
        wdm.cch              = 3;
        
        pls.shape   = 'RC';                       % Shape type
        pls.bw      = 1.0;                        % duty cycle
        pls.ord     = 0;                          % pulse roll-off
        
        snr0dB = zeros(length(pdbm),1);
        
        for k = 1:length(nsc)
            
            sub_signal.nsc       = nsc(k);
            sub_signal.nt        = nt(k);
            sub_signal.nsymb     = 409600/sub_signal.nsc;
            sub_signal.symbrate  = 50/sub_signal.nsc;
            
            signal_prop.symbrate = sub_signal.symbrate*sub_signal.nsc;
            signal_prop.nsymb    = sub_signal.nsymb*sub_signal.nt;
            
            
            
            for i = 1:length(pdbm)
                
                [signals{i},snr0dB(i),ch] = Test_subcarrier(link,sp,signal_prop,sub_signal,amp,pdbm(i),wdm,pls);
                sgs(i).snr0dB = snr0dB;
                sgs(i).sg = 1/sqrt(sub_signal.nsc);
                sgs(i).Pu = pdbm(i);
                
                for j=1:sub_signal.nsc
                    sgs(i).subc(j).tx = sqrt(sub_signal.nsc)*signals{i}.FIELDX_TX(:,j);
                    if (not(sub_signal.nsc == 1))
                        sgs(i).subc(j).rx = signals{i}.SUB_FIELDX(:,j);
                    else
                        sgs(i).subc(j).rx = signals{i}.FIELDX(:,j);
                    end
                end
            end
      
        
        ch_properties       = ch.getProperties;
        ch_properties.Nspan = link.Nspan;
        
        savefile = horzcat('my_sgs_IDA_wdm5_10x100_',int2str(sub_signal.nsc),'sc');
        save(savefile,'sgs','ch_properties','amp','signal_prop','pdbm');
        
        end