
addpath('C:\Users\mene9\Documents\MATLAB\ESSFM_Sim\Test\Safe_Sim\')
% addpath('/home/menelaos/MATLAB/ESSFM_Sim/Test/Safe_Sim/');

symbols      = [2^12];
n_prop_steps = 10;

etasp = [5];
gamma = 1.27e-3;

NS  = [10 10 10];
Nc  = [ 0 0 0];
dBm = (-3:4);

for k = 1:length(etasp)
    for j = 1:length(symbols)
        tic
        for  i = 1:length(NS)
            
            data_3{i} = BER_essfm_vs_ssfm (NS(i),Nc(i),dBm , symbols(j), n_prop_steps,etasp(k),gamma);
            data_1{i} = BER_ESSFM_XY      (NS(i),Nc(i),dBm , symbols(j), n_prop_steps,etasp(k));
            data_2{i} = BER_ESSFM_X       (NS(i),Nc(i),dBm , symbols(j), n_prop_steps,etasp(k));
           
            BER_ESSFM{i} = [data_1{i}(:,2) data_2{i}(:,2) data_3{i}(:,2)];
            
        end
        toc
        
        
        fig = figure(j);
        
        colors  = {'-ob';'-og';'-or';...
            '-+b';'-+g';'-+r';...
            '-*b';'-*g';'-*r'};
        
        lgn = [];
        for i = 1:length(NS)
            
            tmp(i,:)  = ({['ESSFM_{XY}           Ns = ' int2str(NS(i)) ' Nc = ' int2str(Nc(i))],...
                ['ESSFM_X            Ns = ' int2str(NS(i)) ' Nc = ' int2str(Nc(i))],...
                ['Safe ESSFM_X       Ns = ' int2str(NS(i)) ' Nc = ' int2str(Nc(i))]});
            
            lgn           = [lgn tmp(i,:)];
            
        end
        
        
        for i= 1:length(NS)
            p1=semilogy(   dBm , BER_ESSFM{i}(:,1), colors{i}  ,...
                dBm , BER_ESSFM{i}(:,2), colors{i+3},...
                dBm , BER_ESSFM{i}(:,3), colors{i+6});
            hold('on')
        end
        
        t_ESSFM = strcat('propagation of',              {' '},'2^{', int2str(log2(symbols(j))) ,{'}'}, {' '}, ...
            'symbols with'  ,              {' '}, int2str(n_prop_steps)        ,{' '}, ...
            'steps and emission factor of',{' '}, int2str(etasp(k)));
        
        title(t_ESSFM);
        
        grid on;
        ylabel('BER');
        xlabel('Power[dBm]');
        
        legend(lgn)
        
        savefig(fig,strcat('plot/',int2str(log2(symbols(j))),'_',...
            int2str(n_prop_steps),'_',int2str(etasp(k)),'.fig'))
        
        hold('off')
        close(fig);
        
    end
end