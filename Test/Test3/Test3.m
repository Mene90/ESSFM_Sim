
addpath('/home/menelaos/MATLAB/Sim Marco/');

symbols      = 2^14;
n_prop_steps = 10; 

etasp = 4;

NS = [1 8 12];
Nc = [6 4 2];
dBm = (-3:4);

tic 
for i = 1:length(NS)
    
    data_1{i} = BER_ESSFM_Linear_XY(NS(i),Nc(i),dBm , symbols, n_prop_steps,etasp,0);
    data_2{i} = BER_ESSFM_Linear_X (NS(i),Nc(i),dBm , symbols, n_prop_steps,etasp,0);
    data_3{i} = BER_essfm_vs_ssfm  (NS(i),Nc(2),dBm , symbols, n_prop_steps,etasp,0);
    
    BER{i}    = [data_1{i} data_2{i} data_3{i}];
    
end
toc

figure(1);
grid on;
colors  = {'-ob';'-og';'-or';...
           '-+b';'-+g';'-+r';...
           '-*b';'-*g';'-*r'};
       
legends = {strcat('ESSFM_{XY}   Ns =',{' '}, int2str(NS(1)),{' '},'Nc =',{' '}, int2str(Nc(1)));... 
           strcat('ESSFM_X_1    Ns =',{' '}, int2str(NS(1)),{' '},'Nc =',{' '}, int2str(Nc(1)));...
           strcat('ESSFM_X_2    Ns =',{' '}, int2str(NS(1)),{' '},'Nc =',{' '}, int2str(Nc(1)));...
           strcat('ESSFM_{XY}   Ns =',{' '}, int2str(NS(2)),{' '},'Nc =',{' '}, int2str(Nc(2)));...
           strcat('ESSFM_X_1    Ns =',{' '}, int2str(NS(2)),{' '},'Nc =',{' '}, int2str(Nc(2)));...
           strcat('ESSFM_X_2    Ns =',{' '}, int2str(NS(2)),{' '},'Nc =',{' '}, int2str(Nc(2)));...
           strcat('ESSFM_{XY}   Ns =',{' '}, int2str(NS(3)),{' '},'Nc =',{' '}, int2str(Nc(3)));...
           strcat('ESSFM_X_1    Ns =',{' '}, int2str(NS(3)),{' '},'Nc =',{' '}, int2str(Nc(3)));...
           strcat('ESSFM_X_2    Ns =',{' '}, int2str(NS(3)),{' '},'Nc =',{' '}, int2str(Nc(3)))};
       

for i= 1:length(NS)
semilogy(dBm , BER{i}(:,2), colors{i}  ,...
         dBm , BER{i}(:,4), colors{i+3},...
         dBm , BER{i}(:,6), colors{i+6});
hold('on')
end

t = strcat('propagation of',              {' '},'2^{', int2str(log2(symbols)) ,{'}'}, {' '}, ...
           'symbols with'  ,              {' '}, int2str(n_prop_steps)        ,{' '}, ...
           'steps and emission factor of',{' '}, int2str(etasp));
       
title(t)
ylabel('BER');
xlabel('Power[dBm]');
legend(legends{1,1}{1,1},legends{2,1}{1,1},...
       legends{3,1}{1,1},legends{4,1}{1,1},...
       legends{5,1}{1,1},legends{6,1}{1,1},...
       legends{7,1}{1,1},legends{8,1}{1,1},...
       legends{9,1}{1,1});