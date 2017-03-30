% addpath('C:\Users\mene9\Documents\MATLAB\ESSFM_Sim\Test\Safe_Sim\')
% addpath('/home/menelaos/MATLAB/ESSFM_Sim/Test/Safe_Sim/');

<<<<<<< HEAD
symbols      = [2^16];
=======
symbols      = [2^20,2^20,2^20,2^16,2^16];
>>>>>>> 1df47c42ec174171a2ab5fd4d1769672778c6930

forward_steps = [100, 400];
Rs            = [100, 200];
Fn            = [5];
etasp         = [0.5 .*10.^(Fn/10)];
Nspan         = 40;

<<<<<<< HEAD
% N_steps{1}  = [1,2,4,5,8,10,20,40,80,120,160,200]./Nspan;
% N_steps{2}  = [1,2,4,5,8,10,20,40,80,120,160,200]./Nspan;
% N_steps{3}  = [1,2,4,5,8,10,20,40,80,120,160,200]./Nspan;
N_steps{1}  = [1,5,8,10,20,40,80,120,160,200]./Nspan;
N_steps{2}  = [1,5,8,10,20,40,80,120,160,200]./Nspan;

% N_coefficients{1}  = [1,3,17];
% N_coefficients{2}  = [1,3,17];
% N_coefficients{3}  = [1,3,17,23,33,65];
N_coefficients{1}  = [1,3,26,51,101];
N_coefficients{2}  = [1,3,26,51,101];
=======
N_steps{1}  = [1,5,10,20,40,80,120,160,200]./Nspan;
N_steps{2}  = [1,5,10,20,40,80,120,160,200]./Nspan;
N_steps{3}  = [1,5,10,20,40,80,120,160,200]./Nspan;
N_steps{4}  = [1,5,10,20,40,80,120,160,200]./Nspan;
N_steps{5}  = [1,5,10,20,40,80,120,160,200]./Nspan;

N_coefficients{1}  = [1,5,16];
N_coefficients{2}  = [1,5,16];
N_coefficients{3}  = [1,5,17,33];
N_coefficients{4}  = [1,5,33,65];
N_coefficients{5}  = [1,17,33,65,129];
>>>>>>> 1df47c42ec174171a2ab5fd4d1769672778c6930

for k = 1:length(Rs)
    tic
    print = [' '];
    disp(print)
    print = ['Symbrate = ',int2str(Rs(k)),' GBd, Forward steps = ',int2str(forward_steps(k))];
    disp(print);
    print = [' '];
    disp(print);
    
    n_prop_steps = forward_steps(k);
    symbrate     = Rs(k);
    NS           = N_steps{k};
    Nc           = N_coefficients{k};
    
<<<<<<< HEAD
    for  j = 1:length(Nc)
        nc = Nc(j);
        parfor  i = 1:length(NS)
            max_snr(i,j) = ESSFM_MAX_SNR(NS(i),nc,symbols,n_prop_steps,etasp,symbrate,Nspan);
        end
        print = ['ESSFM NS = [',int2str(NS*Nspan),'] Max SNR = [',num2str(max_snr(:,j)'),'] dB ','NC = ',int2str(nc-1)];
        disp(print);
    end
    
=======
    
    disp_comp_max_snr(i) = ones(1,length(NS))*DISPERSION_COMPENSATION_MAXSNR(1,symbols,n_prop_steps,symbrate,etasp,Nspan);
        
    print = ['DISP NS = [',num2str(NS),'] Max SNR = [',num2str(disp_comp_max_snr),'] dB '];
    disp(print);
    
>>>>>>> 1df47c42ec174171a2ab5fd4d1769672778c6930
    
    disp_comp_max_snr = ones(1,length(NS))*DISPERSION_COMPENSATION_MAXSNR(0.025,symbols,n_prop_steps,symbrate,etasp,Nspan);
    
    print = ['DISP NS = [',int2str(NS*Nspan),'] Max SNR = [',num2str(disp_comp_max_snr),'] dB '];
    disp(print);
    
    
%     parfor  i = 1:length(NS)
%         ssfm_max_snr(i) = SSFM_MAX_SNR(NS(i),symbols,n_prop_steps,etasp,symbrate,Nspan);       
%     end
%     
%     print = ['SSFM NS = [',num2str(NS),'] Max SNR = [',num2str(ssfm_max_snr),'] dB '];
%     disp(print);
    
    toc
    
    fig = figure(k);
    
    lgn = [];
    
    colors  = {'-ob';'-og';'-or';'-oc';'-om';'-oc';...
               '-+b';'-+g';'-+r';...
               '-*b';'-*g';'-*r'};
    
    lgn     = [({['DISP']})];
    for i = 1:length(Nc)
        tmp(i,:)  = ({['ESSFM Nc = ' int2str(Nc(i)-1)]});
        lgn       = [lgn tmp(i,:)];
    end
    
    p1 = plot(NS*Nspan, disp_comp_max_snr, '-sk');
    hold('on')
    for i= 1:length(Nc)
        p1=plot(NS*Nspan, max_snr(:,i), colors{i});
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
    xlabel('n� steps');
    
    legend(lgn)
    
    savefig(fig,strcat('plot/S_vs_E_', int2str(n_prop_steps),'_',int2str(symbrate),'_',int2str(log2(symbols)),'.fig'))
    
    hold('off')
    close(fig);
end
