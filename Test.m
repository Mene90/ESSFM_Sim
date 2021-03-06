% nlindex = [1.1e-20,2.5e-20,5e-20];

HPLANCK = 6.62606896e-34;       % Planck's constant [J*s]
CLIGHT = 299792458;             % speed of light [m/s]

n = input('Enter a number: ');

switch n
    case 1
        
%         distribution = input('Enter a distribution: ');
        Nspan            = [12,6,4];  
        link.Nspan       = Nspan(1);
        LL               = [1e5,2e5,3e5];
        link.LL          = LL(1);
        link.attenuation = 0.2;
        link.lambda      = 1550;
        link.sprop       = 5;
        
        sp.bprop    = 1/link.Nspan;
        
        amp.etasp   = 2;
        amp.Fn      = 10*log10(2*amp.etasp);
        
        al =link.attenuation*0.230258509299405*1e-3;
        Gm1=(exp(al*link.LL)-1.d0);
        N0=link.Nspan*Gm1*HPLANCK*CLIGHT/(link.lambda* 1e-9)*amp.etasp;

        signal_prop.nsymb    = 2^24;
        signal_prop.symbrate = 50;
                
%         SNR_dB   = (-10:10:60);
%         pdbm     = 30+10*log10(signal_prop.symbrate*10^9*N0*10.^(SNR_dB*0.1));

        pdbm    =   (-10:5:35);
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
                link.LL    = LL(j);
                sp.bprop    = 1/link.Nspan;
                for i=1:length(pdbm)
                    [signals{i},SNRdB{i},ch] = TestZeroDispHighPower(link,sp,signal_prop,amp,pdbm(i),distribution);
                end
                
                ch_properties       = ch.getProperties;
                ch_properties.Nspan = link.Nspan;
                
                savefile        = strcat('Test_Results/Test1/',distribution,'/',distribution,'_',int2str(link.LL/1000),'X',int2str(link.Nspan),'_NoDisp_DBP_SS');
                
                save(savefile,'signals','SNRdB','ch_properties','amp','signal_prop','pdbm','-v7.3');
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
        signal_prop.nsymb    = 2^20;
        signal_prop.symbrate = 14;     
        
        gpu                  = 0;
        
        for i=1:length(pdbm)
            [signals{i},SNRdB{i},ch] = TestConfronto(link,sp,signal_prop,amp,pdbm(i),'true',gpu);
        end
        
        ch_properties       = ch.getProperties;
        ch_properties.Nspan = link.Nspan;
    
        savefile        = strcat('TestResults/Test2/','G','_',int2str(link.LL/1000),'X',int2str(link.Nspan),'_DBP_Esatta');

        save(savefile,'signals','SNRdB','ch_properties','amp','signal_prop','pdbm','-v7.3');
        
    case 3
        
        distribution     = input('Enter a distribution: ');
        compensation     = [string('ssfm_onlydispcomp'),string('ssfm_onlykercomp'), string('inline')];
        compname         = [string('DispComp'),string('KerrComp'),string('inlineDispComp')];
        Nspan            = (1:5);  
        link.Nspan       = Nspan(1);
        LL               = (0.8e5);
        link.LL          = LL(1);
        link.attenuation = 0.2;
        link.disp        = 16;
        link.n2          = 2.5e-20;
        link.lambda      = 1550;        
        link.sprop       = 25;
        
        sp.bprop    = 25;

        
        amp.etasp   = 1.7;
        amp.type    = 'EDFA';
        amp.Fn      = 10*log10(2*amp.etasp);
        
        al  = link.attenuation*0.230258509299405*1e-3;
        Gm1 = (exp(al*link.LL)-1.d0);
        N0  = link.Nspan*Gm1*HPLANCK*CLIGHT/(link.lambda* 1e-9)*amp.etasp;

        signal_prop.nsymb    = 2^18;
        signal_prop.nt       = 2;
        signal_prop.symbrate = 10;
        
       pls.shape   = 'RRC';                     % Shape type
       pls.bw      = 1;                         % duty cycle
       pls.ord     = 0.2;                       % pulse roll-off
                
%         SNR_dB   = (-10:10:60);
%         pdbm     = 30+10*log10(signal_prop.symbrate*10^9*N0*10.^(SNR_dB*0.1));

        pdbm    =   [2,4,6,8,15];%(-10:1:10);
        P       =   10.^((pdbm-30)*0.1);
        SNR     =   P/signal_prop.symbrate/10^9/N0;
        C       =   log2(1+SNR);
        SNR_dB  =   10*log10(SNR);
        
        gpu     =   1;

        for j = 1:length(pdbm)
            for i=1:length(Nspan)
%                 link.LL          = LL(i);
                link.Nspan       = Nspan(i);
                [signals{i},sig_dbp,SNRdB{i},ch] = TestDispOrKerrComp(link,sp,signal_prop,amp,pdbm(j),distribution,compensation(1),gpu,pls);
            end
            
            ch_properties       = ch{1}.getProperties;
            ch_properties.Nspan = Nspan;
            ch_properties.LL    = LL;
            power.pdbm = pdbm(j);
            
%             savefile        = char(strcat('Test_Results/Test3',distribution,'_',compname(j),'_',int2str(link.LL/1000),'X',int2str(link.Nspan),'_nt_',int2str(signal_prop.nt),'_roll_02'));
            savefile = horzcat('Test_Results/Test3/G_',int2str(link.LL/1000),'X',int2str(link.Nspan),'NoCompensation');
            save(savefile,'signals','SNRdB','ch_properties','amp','signal_prop','pdbm','pls','-v7.3');
%         end
            


            
         end

        
    case 4
        
        link.Nspan       = 10;
        link.LL          = 1e5;
        link.attenuation = 0.2;
        link.lambda      = 1550;
        link.sprop       = 100;
        link.nlindex     = 2.5e-20;
        link.disp        = 17;
        
        sp.bprop    = link.sprop/10;
        
        amp.type    = [string('EDFA'),string('Raman')];
        amp.etasp   = 1;
        amp.Fn      = 10*log10(2*amp.etasp);
        
        al =link.attenuation*0.230258509299405*1e-3;
        Gm1=(exp(al*link.LL)-1.d0);
        N0=link.Nspan*Gm1*HPLANCK*CLIGHT/(link.lambda* 1e-9)*amp.etasp;
        
        pdbm                 = (-15:1:-5);%[-15,-10,-5,-3,-1,0,1,3,5,10];
        signal_prop.nt       = 4;
        signal_prop.nc       = 3;
        signal_prop.nsymb    = 2^20;
        signal_prop.symbrate = 50;  
        
        wdm.cch              = 2;
        
       pls.shape   = 'RRC';                     % Shape type
       pls.bw      = 1;                         % duty cycle
       pls.ord     = 0.2;                       % pulse roll-off
        
        ch_properties.Nspan = link.Nspan;
        gpu = 0;
        
        for j = 1:length(amp.type)
            for i=1:length(pdbm)
                [signals{i},SNRdB{i},ch] = Test_mux(link,sp,signal_prop,amp,pdbm(i),wdm,pls,gpu);
            end
            ch_properties       = ch.getProperties;     
            savefile{j} = horzcat('Test_Results/Test4/',char(amp.type(j)),'/G','_',int2str(link.LL/1000),'X',int2str(link.Nspan),'_WDM_',int2str(signal_prop.nc),'_',char(amp.type(j)),'_nt_',int2str(signal_prop.nt),'_roll_02');
            save(savefile{j},'signals','SNRdB','ch_properties','amp','signal_prop','pdbm','pls','-v7.3');
% =======
%             savefile(j) = char(strcat('Test_Results/Test4/',amp.type(j),'/G','_',int2str(link.LL/1000),'X',int2str(link.Nspan),'_WDM_',int2str(signal_prop.nc),'_',amp.type(j),'_nt_',int2str(signal_prop.nt),'roll_01'));
%             save(savefile(j),'signals','SNRdB','ch_properties','amp','signal_prop','pdbm','pls','-v7.3');
% >>>>>>> 54b6c003a0546510bdcf15474094d40eb801016f
        end
        
         
%       
%         
%         for i =1:length(pdbm)
%             irate(i) = AIR(signals{i},'Gaussian');
%         end
%         
%         ylabel('IR [bit/symbol]')
%         xlabel('launch power per channel [dBm]')
%         grid on
%         plot(pdbm,irate,'-ob');
        
%         savefile = strcat('Test_Results/Test4/',amp.type,'/G','_',int2str(link.LL/1000),'X',int2str(link.Nspan),'_WDM_',int2str(signal_prop.nc),'_nt_',int2str(signal_prop.nt));

%         save(savefile,'signals','SNRdB','ch_properties','amp','signal_prop','pdbm');
        
end