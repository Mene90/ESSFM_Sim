% addpath('C:\Users\mene9\Documents\MATLAB\ESSFM_Sim\Test\Safe_Sim\')
addpath('/home/menelaos/MATLAB/ESSFM_Sim/Test/Safe_Sim/');

symbols      = [2^14];
n_prop_steps = 20;

symbrate = 32;
etasp = [1.5];
Nspan = 40;

NS  = [1,2,4,6];
Nc  = [1,8];

for  j = 1:length(Nc)
    for  i = 1:length(NS)
        
        max_snr(i,j) = ESSFM_MAX_SNR(NS(i),Nc(j),symbols,n_prop_steps,symbrate,etasp,Nspan);
        
    end
end
fig = figure(j);

lgn = [];

colors  = {'-ob';'-og';'-or';...
           '-+b';'-+g';'-+r';...
           '-*b';'-*g';'-*r'};
for i = 1:length(Nc)
    tmp(i,:)  = ({['ESSFM Nc = ' int2str(Nc(i))]});
    lgn       = [lgn tmp(i,:)];
end

for i= 1:length(Nc)
    p1=plot(NS, max_snr(:,i), colors{i});
    hold('on')
end

title('SSFM vs ESSFM');

grid on;
ylabel('SNR [dB]');
xlabel('n° steps');

legend(lgn)

savefig(fig,strcat('plot/', int2str(n_prop_steps),'_',int2str(Nc(j)),'.fig'))

hold('off')
close(fig);

