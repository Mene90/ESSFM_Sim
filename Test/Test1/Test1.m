
addpath('/home/menelaos/MATLAB/Sim Marco/');

symbols      = 2^12;
n_prop_steps = 100; 

NS = [2 4 12];
Nc = [6 4 2];
dBm = (-3:4);

etasp = 4;

tic 
for i = 1:length(NS)
    data_1{i} = BER_ESSFM_vs_SSFM(NS(i),Nc(i),dBm , symbols, n_prop_steps,etasp);
    data_2{i} = BER_essfm_vs_ssfm(NS(i),Nc(i),dBm , symbols, n_prop_steps,etasp);
    BER{i}   = [data_1{i} data_2{i}];
end
toc

figure(1);
grid on;
colors  = {'-ob';'-og';'-or';'-+b';'-+g';'-+r'};
legends = {strcat('ESSFM_1  Ns =',{' '}, int2str(NS(1)),{' '},'Nc =',{' '}, int2str(Nc(1)));... 
           strcat('ESSFM_2  Ns =',{' '}, int2str(NS(1)),{' '},'Nc =',{' '}, int2str(Nc(1)));...
           strcat('ESSFM_1  Ns =',{' '}, int2str(NS(2)),{' '},'Nc =',{' '}, int2str(Nc(2)));...
           strcat('ESSFM_2  Ns =',{' '}, int2str(NS(2)),{' '},'Nc =',{' '}, int2str(Nc(2)));...
           strcat('ESSFM_1  Ns =',{' '}, int2str(NS(3)),{' '},'Nc =',{' '}, int2str(Nc(3)));...
           strcat('ESSFM_2  Ns =',{' '}, int2str(NS(3)),{' '},'Nc =',{' '}, int2str(Nc(3)))};

for i= 1:length(NS)
semilogy(dBm , BER{i}(:,2), colors{i}, dBm , BER{i}(:,4), colors{i+3});
hold('on')
end
t = strcat('Propagaion of',{' '}, int2str(symbols)     ,{' '}, ...
           'symbols with' ,{' '}, int2str(n_prop_steps),{' '}, 'steps');
title(t)
ylabel('BER');
xlabel('Power[dBm]');
legend(legends{1,1}{1,1},legends{2,1}{1,1},...
       legends{3,1}{1,1},legends{4,1}{1,1},...
       legends{5,1}{1,1},legends{6,1}{1,1});



