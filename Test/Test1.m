
addpath('/home/menelaos/MATLAB/Sim Marco/');

% fscanformat = '%f %f %f %f';
% fprintformat = '%f %f %f %f\r\n';

% if(exist('test1_result', 'file'))
%     fileID = fopen('test1_result','r');
% else
%     fileID = fopen('test1_result','w+');
% end

symbols      = 2^16;
n_prop_steps = 100; 

NS = [1 2 8];
Nc = [1 2 2];
dBm = (-3:4);
for i = 1:length(NS)
    data_1 = BER_ESSFM_vs_SSFM(NS(i),Nc(i),dBm , symbols, n_prop_steps);
    data_2 = BER_essfm_vs_ssfm(NS(i),Nc(i),dBm , symbols, n_prop_steps);
    BER{i}   = [data_1 data_2];
%     fprintf(fileID,fprintformat , BER{i});
end


figure(1);
grid on;
colors  = {'-ob';'-og';'-or';'-+b';'-+g';'-+r'};
legends = {'ESSFM_1 Ns = 1 Nc = 1' ; 'ESSFM_2 Ns = 1 Nc = 1';...
           'ESSFM_1 Ns = 2 Nc = 2' ; 'ESSFM_2 Ns = 2 Nc = 2';...
           'ESSFM_1 Ns = 8 Nc = 2' ; 'ESSFM_2 Ns = 8 Nc = 2'};

for i= 1:length(NS)
semilogy(dBm , BER{i}(:,2), colors{i}, dBm , BER{i}(:,4), colors{i+3});
hold('on')
end

ylabel('BER');
xlabel('Power[dBm]');
legend(legends{1},legends{2},legends{3},legends{4},legends{5},legends{6});



