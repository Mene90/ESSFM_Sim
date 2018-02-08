%         clc;
        clear all;

        HPLANCK = 6.62606896e-34;       % Planck's constant [J*s]
        CLIGHT = 299792458;             % speed of light [m/s]
      
        
        Nspan            = [100];
        link.LL          = 0.6e5;
        link.attenuation = 0.2;
        link.lambda      = 1550;
        link.sprop       = 1000;
        link.nlindex     = 2.5e-20;
        link.disp        = 17;
        
        sp.bprop    = 100;
        
        amp.type    = 'EDFA';
        amp.etasp   = 1.6;
        amp.Fn      = 10*log10(2*amp.etasp);
        
%         al =link.attenuation*0.230258509299405*1e-3;
%         Gm1=(exp(al*link.LL)-1.d0);
%         N0=link.Nspan*Gm1*HPLANCK*CLIGHT/(link.lambda* 1e-9)*amp.etasp;
        
        nsc                  = [1,2,4];
        nt                   = [1,2,4];       
        
        pdbm                 = (-6:1:0);
        
        signal_prop.nt       = 8;
        signal_prop.nc       = 5;  
        
        pol                  = 1;
   
        
        wdm.cch              = 3;
        
        pls.shape   = 'RC';                       % Shape type
        pls.bw      = 1.0;                        % duty cycle
        pls.ord     = 0;                          % pulse roll-off
        
        SymbolsXSubcarrier = 204800;
        
        gpu = 1;
       
        for m = 1:length(Nspan)
            disp(Nspan(m));
            tic
            for k = 1:length(nsc)
                
                
                link.Nspan = Nspan(m);
                
                snr0dB = zeros(length(pdbm),1);
                
                sub_signal.nsc       = nsc(k);
                sub_signal.nt        = nt(k);
                sub_signal.nsymb     = SymbolsXSubcarrier;
                sub_signal.symbrate  = 50/sub_signal.nsc;
                
                signal_prop.symbrate = sub_signal.symbrate*sub_signal.nsc;
                signal_prop.nsymb    = sub_signal.nsymb*sub_signal.nt;
                
                parfor i = 1:length(pdbm)
                    
                    [signals{i},snr0dB(i),ch{i}] = Test_subcarrier(link,sp,signal_prop,sub_signal,amp,pdbm(i),wdm,pls,pol,gpu);
                    
                    sgs(i).sg = 1/sqrt(sub_signal.nsc);
                    sgs(i).Pu = pdbm(i);
                    
                    
                    for j=1:sub_signal.nsc
                        
                        
                        n_column = 200;
                        rowsize = SymbolsXSubcarrier/n_column;
                        
                        if ( pol == 1 )
                            sgs(i).subc(j).tx = reshape(sqrt(sub_signal.nsc)*signals{i}.FIELDX_TX(:,j),rowsize,n_column,1);
                            if (not(sub_signal.nsc == 1))
                                sgs(i).subc(j).rx = reshape(signals{i}.SUB_FIELDX(:,j),rowsize,n_column,1);
                            else
                                sgs(i).subc(j).rx = reshape(signals{i}.FIELDX(:,j),rowsize,n_column,1);
                            end
                        else
                            sgs(i).subc(j).tx = reshape(sqrt(sub_signal.nsc)*[signals{i}.FIELDX_TX(:,j);signals{i}.FIELDY_TX(:,j)],rowsize,n_column,2);
                            if (not(sub_signal.nsc == 1))
                                sgs(i).subc(j).rx = reshape([signals{i}.SUB_FIELDX(:,j);signals{i}.SUB_FIELDY(:,j)],rowsize,n_column,2);
                            else
                                sgs(i).subc(j).rx = reshape([signals{i}.FIELDX(:,j);signals{i}.FIELDY(:,j)],rowsize,n_column,2);
                            end
                        end
                    end
                    
                    
                end
               
                for i = 1:length(pdbm)
                    sgs(i).snr0dB = zeros(length(pdbm),1);
                    sgs(i).snr0dB(1:i) = snr0dB(1:i);
                end
                
                ch_properties       = ch{1}.getProperties;
                ch_properties.Nspan = link.Nspan;
                if ( pol == 2 )
                    savefile = horzcat('my_sgs_dp_LA_wdm5_',int2str(Nspan(m)),'x60_',int2str(sub_signal.nsc),'sc');
                else
                    savefile = horzcat('my_sgs_LA_wdm5_',int2str(Nspan(m)),'x60_',int2str(sub_signal.nsc),'sc');
                end
                save(savefile,'sgs','ch_properties','amp','signal_prop','pdbm','-v7.3');
                clear signals ;
            end
            t = toc;
            disp(horzcat('Total time: ',num2str(round(t/60/60,2)),'hours'));
            
        end