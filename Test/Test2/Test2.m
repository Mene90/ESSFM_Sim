
% addpath('C:\Users\mene9\Documents\MATLAB\ESSFM_Sim\Test\Safe_Sim\')
% addpath('/home/menelaos/MATLAB/ESSFM_Sim/Test/Safe_Sim/');

symbols      = [2^18];
n_prop_steps = 10;

Fn       = [7];
etasp    = [0.5 .*10.^(Fn/10)];
gamma = 1.2e-3;

NS  = [0.025,0.025,0.025];
Nc  = [1  5 45];
dBm = (-2:6);

for j = 1:length(symbols)
    
    for i = 1:length(NS)
        
        data_1{i} = BER_ESSFM_XY(NS(i),Nc(i),dBm , symbols(j), n_prop_steps,etasp); 
        data_3{i} = BER_essfm_vs_ssfm (NS(i),Nc(i),dBm , symbols(j), n_prop_steps,etasp,gamma); 
        data_2{i} = BER_ESSFM_X (NS(i),Nc(i),dBm , symbols(j), n_prop_steps,etasp);
               
        BER{i}    = [data_1{i} data_2{i} data_3{i}];
        
    end
    
    
    
    fig = figure(j);
    colors  = {'-ob';'-og';'-or';...
        '-+b';'-+g';'-+r';...
        '-*b';'-*g';'-*r'};
    
    legends = { strcat('ESSFM_{XY}   Ns =',{' '}, int2str(NS(1)*40),{' '},'Nc =',{' '}, int2str(Nc(1)-1));...
        strcat('ESSFM_X_1            Ns =',{' '}, int2str(NS(1)*40),{' '},'Nc =',{' '}, int2str(Nc(1)-1));...
        strcat('Safe_ESSFM_X         Ns =',{' '}, int2str(NS(1)*40),{' '},'Nc =',{' '}, int2str(Nc(1)-1));...
        strcat('ESSFM_{XY}           Ns =',{' '}, int2str(NS(2)*40),{' '},'Nc =',{' '}, int2str(Nc(2)-1));...
        strcat('ESSFM_X_1            Ns =',{' '}, int2str(NS(2)*40),{' '},'Nc =',{' '}, int2str(Nc(2)-1));...
        strcat('Safe_ESSFM_X         Ns =',{' '}, int2str(NS(2)*40),{' '},'Nc =',{' '}, int2str(Nc(2)-1));...
        strcat('ESSFM_{XY}           Ns =',{' '}, int2str(NS(3)*40),{' '},'Nc =',{' '}, int2str(Nc(3)-1));...
        strcat('ESSFM_X_1            Ns =',{' '}, int2str(NS(3)*40),{' '},'Nc =',{' '}, int2str(Nc(3)-1));...
        strcat('Safe_ESSFM_X         Ns =',{' '}, int2str(NS(3)*40),{' '},'Nc =',{' '}, int2str(Nc(3)-1))};
    
    
    for i= 1:length(NS)
        semilogy(   dBm , BER{i}(:,2), colors{i}  ,...
            dBm , BER{i}(:,4), colors{i+3},...
            dBm , BER{i}(:,6), colors{i+6});
            
            
        hold('on')
    end
    
    t = strcat('propagation of',              {' '},'2^{', int2str(log2(symbols(j))) ,{'}'}, {' '}, ...
        'symbols with'  ,              {' '}, int2str(n_prop_steps)        ,{' '}, ...
        'steps and emission factor of',{' '}, int2str(etasp));
    
    title(t)
    
    grid on;
    ylabel('BER');
    xlabel('Power[dBm]');
    legend(legends{1,1}{1,1},legends{2,1}{1,1},...
        legends{3,1}{1,1},legends{4,1}{1,1},...
        legends{5,1}{1,1},legends{6,1}{1,1},...
        legends{7,1}{1,1},legends{8,1}{1,1},...
        legends{9,1}{1,1});
    
    savefig(fig,strcat('plot/',int2str(log2(symbols(j))),'_',...
                      int2str(n_prop_steps),'_',int2str(etasp),'.fig'))
    hold('off')
    close(fig);
   
end
