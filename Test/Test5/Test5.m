% addpath('C:\Users\mene9\Documents\MATLAB\ESSFM_Sim\Test\Safe_Sim\')
% addpath('/home/menelaos/MATLAB/ESSFM_Sim/Test/Safe_Sim/');

symbols      = [2^16];
n_prop_steps = 400;

symbrate = 200;
Fn       = [5];
etasp    = [0.5 .*10.^(Fn/10)];
Nspan    = 40;

NS  = [0.025,0.25,0.5,1,2,4,8,25,50,100];
Nc  = [1,2,3,5,17,33];

tic
for  j = 1:length(Nc)
    nc = Nc(j);
    for  i = 1:length(NS)
        max_snr(i,j) = ESSFM_MAX_SNR(NS(i),nc,symbols,n_prop_steps,symbrate,etasp,Nspan);
    end
    print = ['ESSFM NS = [',int2str(NS),'] Max SNR = [',int2str(max_snr(:,j)'),'] dB ','NC = ',int2str(nc)];
    disp(print);
end
toc

tic
for  i = 1:length(NS)
    disp_comp_max_snr(i) = DISP_COMP_MAX_SNR(NS(i),symbols,n_prop_steps,symbrate,etasp,Nspan);
end

print = ['DISP NS = [',int2str(NS),'] Max SNR = [',int2str(disp_comp_max_snr),'] dB '];
disp(print);
toc

tic
for  i = 1:length(NS)
        ssfm_max_snr(i) = SSFM_MAX_SNR(NS(i),symbols,n_prop_steps,symbrate,etasp,Nspan);        
end
print = ['SSFM NS = [',int2str(NS),'] Max SNR = [',int2str(ssfm_max_snr),'] dB '];
disp(print);
toc


fig = figure(j);

lgn = [];

colors  = {'-ob';'-og';'-or';'-oc';'-om';'-oc';...
           '-+b';'-+g';'-+r';...
           '-*b';'-*g';'-*r'};
       
lgn     = [({['DISP']}) ({['SSFM']})]
for i = 1:length(Nc)
    tmp(i,:)  = ({['ESSFM Nc = ' int2str(Nc(i))]});
    lgn       = [lgn tmp(i,:)];
end

p1 = plot(NS, disp_comp_max_snr, '-ok', NS, ssfm_max_snr, '-*k');
hold('on')
for i= 1:length(Nc)
    p1=plot(NS, max_snr(:,i), colors{i});
    hold('on')
end
t = strcat('SSFM vs ESSFM N_{step} =',{' '},int2str(n_prop_steps),{' '},...
           'F_n = '    ,{' '}, int2str(Fn)    ,'dB',{' '},...
           'N_{Span}= ',{' '}, int2str(Nspan) ,     {' '},...
           'R = '      ,{' '}, int2str(symbrate),   {' '},...
           'N_{symbols}',{' '}, '2^{', int2str(log2(symbols)) ,{'}'});
       
title(t);

grid on;
ylabel('SNR [dB]');
xlabel('n° steps');

legend(lgn)

savefig(fig,strcat('plot/S_vs_E_', int2str(n_prop_steps),'_',int2str(symbrate),'_',int2str(log2(symbols)),'.fig'))

hold('off')
close(fig);

