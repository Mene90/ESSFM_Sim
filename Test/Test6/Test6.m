% addpath('C:\Users\mene9\Documents\MATLAB\ESSFM_Sim\Test\Safe_Sim\')
% addpath('/home/menelaos/MATLAB/ESSFM_Sim/Test/Safe_Sim/');

symbols      = [2^16];
n_prop_steps = 50;

symbrate = 50;
Fn       = [5];
etasp    = [0.5 .*10.^(Fn/10)];
Nspan    = [10:10:80];

NS  = [1,2,4,5];
Nc  = [13];


for  k = 1:length(NS)
    tic
    
    for i = 1:length(Nspan)
        ssfm_max_snr(k,i) = ESSFM_MAX_SNR(NS(k)/Nspan(i),1,symbols,n_prop_steps,etasp,symbrate,Nspan(i));
    end
 
    print = ['SSFM  Nspan = [',int2str(Nspan),'] Max SNR = [',int2str(ssfm_max_snr(k,:)),'] dB NS = ', int2str(NS(k))];
    disp(print);
      
    for  i = 1:length(Nspan)
        max_snr(k,i) = ESSFM_MAX_SNR(NS(k)/Nspan(i),Nc,symbols,n_prop_steps,etasp,symbrate,Nspan(i));
    end
    
    print = ['ESSFM Nspan = [',int2str(Nspan),'] Max SNR = [',int2str(max_snr(k,:)),'] dB ','NS = ',int2str(NS(k))];
    disp(print);
    toc
    
end

parfor  i = 1:length(Nspan)
    disp_comp_max_snr(1,i) = DISPERSION_COMPENSATION_MAXSNR(NS(1)/Nspan(i),symbols,n_prop_steps,symbrate,etasp,Nspan(i));
end

print = ['DISP NS = [',int2str(Nspan),'] Max SNR = [',int2str(disp_comp_max_snr(1,:)),'] dB '];
disp(print);

fig = figure(1);

lgn = [];

colors  = {'-or' '-+r' '-*r' '-sr' '-dr';...
           '-og' '-+g' '-*g' '-sg' '-dg'};
      
       
for i = 1:length(NS)
    tmp(i,:)  = ({['SSFM N_{steps} = ' int2str(NS(i))]});
    lgn       = [lgn tmp(i,:)];
end 

for i = 1:length(NS)
    tmp(i,:)  = ({['ESSFM N_{steps} = ' int2str(NS(i))]});
    lgn       = [lgn tmp(i,:)];
end

for i = 1:length(NS)
    tmp(i,:)  = ({['SSFM DISP. COMP = ' int2str(NS(i))]});
    lgn       = [lgn tmp(i,:)];
end 


for i= 1:length(NS)
    plot(Nspan, ssfm_max_snr(i,:), colors{1,i});
    hold('on')
end

for i= 1:length(NS)
    plot(Nspan, max_snr(i,:), colors{2,i});
    hold('on')
end


p1=plot(Nspan, disp_comp_max_snr(1,1:8));
hold('on')


t = strcat('SNR vs N_{span} N_{step} =',{' '},int2str(n_prop_steps),{' '},...
           'F_n = '     ,{' '}, int2str(Fn)    ,'dB',{' '},...
           'R_s = '     ,{' '}, int2str(symbrate),   {' '},...
           'N_{symbols}',{' '}, '2^{', int2str(log2(symbols)) ,{'}'});
       
title(t);

grid on;
ylabel('SNR [dB]');
xlabel('Distance [n� Span]');

legend(lgn)

savefig(fig,strcat('plot/SNR_Nspan_', int2str(n_prop_steps),'_',int2str(symbrate),'_',int2str(log2(symbols)),'.fig'))

hold('off')
close(fig);