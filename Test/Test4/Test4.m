% addpath('C:\Users\mene9\Documents\MATLAB\ESSFM_Sim\Test\Safe_Sim\')
addpath('/home/menelaos/MATLAB/ESSFM_Sim/Test/Safe_Sim/');

symbols      = [2^16];
n_prop_steps = 10;

etasp = [4];
gamma = 1.27e-3;

NS  = [1];
Nc  = [ 0];


for j = 1:length(symbols)
    tic
    for  i = 1:length(NS)
        
        max = MIN_BER_ESSFM_XY(NS(i),Nc(i),symbols(j),n_prop_steps,etasp);
        
    end
    toc
%     
%     
%     
%             fig = figure(j);
%     
%             colors  = {'-ob';'-og';'-or';...
%                 '-+b';'-+g';'-+r';...
%                 '-*b';'-*g';'-*r'};
%     
%             lgn = [];
%             for i = 1:length(NS)
%     
%                 tmp(i,:)  = ({['ESSFM_{XY}           Ns = ' int2str(NS(i)) ' Nc = ' int2str(Nc(i))],...
%                     ['ESSFM_X            Ns = ' int2str(NS(i)) ' Nc = ' int2str(Nc(i))],...
%                     ['Safe ESSFM_X       Ns = ' int2str(NS(i)) ' Nc = ' int2str(Nc(i))]});
%     
%                 lgn           = [lgn tmp(i,:)];
%     
%             end
%     
%     
%             for i= 1:length(NS)
%                 p1=semilogy(   dBm , BER_ESSFM{i}(:,1), colors{i}  ,...
%                     dBm , BER_ESSFM{i}(:,2), colors{i+3},...
%                     dBm , BER_ESSFM{i}(:,3), colors{i+6});
%                 hold('on')
%             end
%     
%             t_ESSFM = strcat('propagation of',              {' '},'2^{', int2str(log2(symbols(j))) ,{'}'}, {' '}, ...
%                 'symbols with'  ,              {' '}, int2str(n_prop_steps)        ,{' '}, ...
%                 'steps and emission factor of',{' '}, int2str(etasp(k)));
%     
%             title(t_ESSFM);
%     
%             grid on;
%             ylabel('BER');
%             xlabel('Power[dBm]');
%     
%             legend(lgn)
%     
%             savefig(fig,strcat('plot/',int2str(log2(symbols(j))),'_',...
%                 int2str(n_prop_steps),'_',int2str(etasp(k)),'.fig'))
%     
%             hold('off')
%             close(fig);
    
end
