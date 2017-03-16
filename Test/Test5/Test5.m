% addpath('C:\Users\mene9\Documents\MATLAB\ESSFM_Sim\Test\Safe_Sim\')
% addpath('/home/menelaos/MATLAB/ESSFM_Sim/Test/Safe_Sim/');

symbols      = [2^14];
n_prop_steps = 10;

symbrate = 32;
Fn       = [5];
etasp    = [0.5 .*10.^(Fn/10)];
Nspan    = 40;
% gamma    = 1.27e-3;

NS  = [1,2,4,5,10,20,40,80]./Nspan;
Nc  = [1,5];

tic
for  j = 1:length(Nc)
    nc = Nc(j);
    parfor  i = 1:length(NS)
        max_snr(i,j) = ESSFM_MAX_SNR(NS(i),nc,symbols,n_prop_steps,etasp,symbrate,Nspan);
    end
    print = ['ESSFM NS = [',int2str(NS.*Nspan),'] Max SNR = [',num2str(max_snr(:,j)'),'] dB ','NC = ',int2str(nc-1)];
    disp(print);
end
toc

tic
disp_comp_max_snr = ones(1,length(NS))*DISPERSION_COMPENSATION_MAXSNR(1,symbols,n_prop_steps,symbrate,etasp,Nspan);
print = ['DISP NS = [',int2str(NS.*Nspan),'] Max SNR = [',int2str(disp_comp_max_snr),'] dB '];
disp(print);
toc

% tic
% parfor  i = 1:length(NS)
%         ssfm_max_snr(i) = SSFM_MAX_SNR(NS(i),symbols,n_prop_steps,symbrate,etasp,Nspan);        
% end
% print = ['SSFM NS = [',int2str(NS.*Nspan),'] Max SNR = [',int2str(ssfm_max_snr),'] dB '];
% disp(print);
% toc


fig = figure(j);

lgn = [];

colors  = {'-ob';'-og';'-or';'-oc';'-om';'-oc';...
           '-+b';'-+g';'-+r';...
           '-*b';'-*g';'-*r'};
       

% lgn     = [({['Disp. Comp.']})];
for i = 1:length (Nc)
    tmp(i,:)  = ({['ESSFM Nc = ' int2str(Nc(i)-1)]});
    lgn       = [lgn tmp(i,:)];
end

p1 = plot(NS.*Nspan, disp_comp_max_snr, '--k');
hold('on')
for i= 1:length(Nc)
    p1=plot(NS.*Nspan, max_snr(1:10,i), colors{i});
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

savefig(fig,strcat('plot/ESSFM_', int2str(n_prop_steps),'_',int2str(symbrate),'_',int2str(log2(symbols)),'.fig'))

hold('off')
close(fig);

