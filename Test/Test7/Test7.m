% addpath('C:\Users\mene9\Documents\MATLAB\ESSFM_Sim\Test\Safe_Sim\')
% addpath('/home/menelaos/MATLAB/ESSFM_Sim/Test/Safe_Sim/');

symbols      = [2^18];

forward_steps = [10, 10,  50, 200, 800];
Rs            = [10, 25,  50, 100, 200];
Fn            = [5];
etasp         = [0.5 .*10.^(Fn/10)];
Nspan         = 40;

N_steps{1}  = [1,2,4,5,8,10,20,40,80,120,160,200]./Nspan;
N_steps{2}  = [1,2,4,5,8,10,20,40,80,120,160,200]./Nspan;
N_steps{3}  = [1,2,4,5,8,10,20,40,80,120,160,200]./Nspan;
N_steps{4}  = [1,2,4,5,8,10,20,40,80,120,160,200]./Nspan;
N_steps{5}  = [1,2,4,5,8,10,20,40,80,120,160,200]./Nspan;

N_coefficients{1}  = [1,3,5];
N_coefficients{2}  = [1,3,5];
N_coefficients{3}  = [1,3,5,9,13,17];
N_coefficients{4}  = [1,3,5,9,13,17];
N_coefficients{5}  = [1,3,5,9,17,33];

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
    
    parfor  i = 1:length(NS)
        disp_comp_max_snr(i) = DISPERSION_COMPENSATION_MAXSNR(NS(i),symbols,n_prop_steps,symbrate,etasp,Nspan);
    end
    
    print = ['DISP NS = [',num2str(NS),'] Max SNR = [',num2str(disp_comp_max_snr),'] dB '];
    disp(print);
    
    
    parfor  i = 1:length(NS)
        ssfm_max_snr(i) = SSFM_MAX_SNR(NS(i),symbols,n_prop_steps,symbrate,etasp,Nspan);       
    end
    
    print = ['SSFM NS = [',num2str(NS),'] Max SNR = [',num2str(ssfm_max_snr),'] dB '];
    disp(print);
    
    
    
    for  j = 1:length(Nc)
        nc = Nc(j);
        parfor  i = 1:length(NS)
            max_snr(i,j) = ESSFM_MAX_SNR(NS(i),nc,symbols,n_prop_steps,symbrate,etasp,Nspan);
        end
        print = ['ESSFM NS = [',num2str(NS),'] Max SNR = [',num2str(max_snr(:,j)'),'] dB ','NC = ',int2str(nc-1)];
        disp(print);
    end
    toc
    
    fig = figure(k);
    
    lgn = [];
    
    colors  = {'-ob';'-og';'-or';'-oc';'-om';'-oc';...
               '-+b';'-+g';'-+r';...
               '-*b';'-*g';'-*r'};
    
    lgn     = [({['DISP']}) ({['SSFM']})];
    for i = 1:length(Nc)
        tmp(i,:)  = ({['ESSFM Nc = ' int2str(Nc(i)-1)]});
        lgn       = [lgn tmp(i,:)];
    end
    
    p1 = plot(NS, disp_comp_max_snr, '-sk', NS, ssfm_max_snr, '-*k');
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
end
