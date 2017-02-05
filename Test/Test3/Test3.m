% addpath('C:\Users\mene9\Documents\MATLAB\ESSFM_Sim\Test\Safe_Sim\')
addpath('/home/menelaos/MATLAB/ESSFM_Sim/Test/Safe_Sim/');

symbols      = [2^12];
n_prop_steps = 10;

etasp = 0;

NS  = [1 8 12];
Nc  = [1 1  1];
dBm = (-3:4);

for j = 1:length(etasp)
    for i = 1:length(NS)
        
        data_1{i} = BER_ESSFM_Linear_XY(NS(i),Nc(i),dBm , symbols(j), n_prop_steps,etasp,0);
        data_2{i} = BER_ESSFM_Linear_X (NS(i),Nc(i),dBm , symbols(j), n_prop_steps,etasp,0);        
        BER{i}    = [data_1{i} data_2{i}];
        
    end
    
    
    
    fig = figure(j);
    colors  = {'-ob';'-og';'-or';...
        '-+b';'-+g';'-+r';...
        '-*b';'-*g';'-*r'};
    
    legends = { strcat('ESSFM_{XY}   Ns =',{' '}, int2str(NS(1)),{' '},'Nc =',{' '}, int2str(Nc(1)));...
            strcat('ESSFM_X      Ns =',{' '}, int2str(NS(1)),{' '},'Nc =',{' '}, int2str(Nc(1)));...
            strcat('ESSFM_{XY}   Ns =',{' '}, int2str(NS(2)),{' '},'Nc =',{' '}, int2str(Nc(2)));...
            strcat('ESSFM_X      Ns =',{' '}, int2str(NS(2)),{' '},'Nc =',{' '}, int2str(Nc(2)));...
            strcat('ESSFM_{XY}   Ns =',{' '}, int2str(NS(3)),{' '},'Nc =',{' '}, int2str(Nc(3)));...
            strcat('ESSFM_X      Ns =',{' '}, int2str(NS(3)),{' '},'Nc =',{' '}, int2str(Nc(3)))};
    
    
    for i= 1:length(NS)
        plot(   dBm , BER{i}(:,2), colors{i}  ,...
            dBm , BER{i}(:,4), colors{i+3});
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
        legends{5,1}{1,1},legends{6,1}{1,1});
    
    savefig(fig,strcat('plot/',int2str(log2(symbols(j))),'_',...
                      int2str(n_prop_steps),'_',int2str(etasp),'.fig'))
    hold('off')
    close(fig);
   
end