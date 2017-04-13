% addpath('C:\Users\mene9\Documents\MATLAB\ESSFM_Sim\Test\Safe_Sim\')
% addpath('/home/menelaos/MATLAB/ESSFM_Sim/Test/Safe_Sim/');

symbols      = [2^16];
n_prop_steps = 100;

symbrate = 50;
Fn       = [5];
etasp    = [0.5 .*10.^(Fn/10)];
Nspan    = [10:10:80];

NS  = [1,2,3,4,5,6,7,8];
Nc  = [65];


% for  k = 1:length(NS)
    tic
    
%     parfor  i = 1:length(Nspan)
%         max_snr(1,i) = ESSFM_MAX_SNR(NS(i)/Nspan(i),Nc,symbols,n_prop_steps,etasp,symbrate,Nspan(i),[-8 8]);
%     end
    
    print = ['ESSFM Nspan = [',int2str(Nspan),'] Max SNR = [',num2str(max_snr(1,:)),'] dB ','NS = [',int2str(NS),']'];
    disp(print);
    
    parfor i = 1:length(Nspan)
        ssfm_max_snr(1,i) = ESSFM_MAX_SNR(NS(i)/Nspan(i),1,symbols,n_prop_steps,etasp,symbrate,Nspan(i),[-8 8]);
    end
 
    print = ['SSFM  Nspan = [',int2str(Nspan),'] Max SNR = [',int2str(ssfm_max_snr(1,:)),'] dB NS = [',int2str(NS),']'];
    disp(print);
      
    toc
    
% end

parfor  i = 1:length(Nspan)
    disp_comp_max_snr(1,i) = DISPERSION_COMPENSATION_MAXSNR(NS(i)/Nspan(i),symbols,n_prop_steps,symbrate,etasp,Nspan(i));
end

print = ['DISP NS = [',int2str(Nspan),'] Max SNR = [',int2str(disp_comp_max_snr(1,:)),'] dB '];
disp(print);

fig = figure(1);

lgn = [];

colors  = {'-or' '-+r' '-*r' '-sr' '-dr';...
           '-og' '-+g' '-*g' '-sg' '-dg'};
      
       

    tmp(1,:)  = ({['SSFM N_{steps} = ' int2str(NS(1)*100*10) 'km']});
    lgn       = [lgn tmp(1,:)];


    tmp(1,:)  = ({['ESSFM N_{steps} = ' int2str(NS(1)*100*10) 'km']});
    lgn       = [lgn tmp(1,:)];



    tmp(1,:)  = ({['DISP. COMP']});
    lgn       = [lgn tmp(1,:)];



% for i= 1:length(NS)
    plot(Nspan, ssfm_max_snr(1,:), colors{1,1});
    hold('on')
% end

% for i= 1:length(NS)
    plot(Nspan, max_snr(1,:), colors{2,1});
    hold('on')
% end


p1=plot(Nspan, disp_comp_max_snr(1,:));
hold('on')


t = strcat('SNR vs N_{span} N_{step} =',{' '},int2str(n_prop_steps),{' '},...
           'F_n = '     ,{' '}, int2str(Fn)    ,'dB',{' '},...
           'R_s = '     ,{' '}, int2str(symbrate),   {' '},...
           'N_{symbols}',{' '}, '2^{', int2str(log2(symbols)) ,{'}'});
       
title(t);

grid on;
ylabel('SNR [dB]');
xlabel('Distance [n° Span]');

legend(lgn)

savefig(fig,strcat('plot/SNR_Nspan_', int2str(n_prop_steps),'_',int2str(symbrate),'_',int2str(log2(symbols)),'_1000kmXstep.fig'))

hold('off')
close(fig);