symbols      = [2^16];
n_prop_steps = 50;

Fn       = [7];
etasp    = [0.5 .*10.^(Fn/10)];
gamma = 1.2e-3;
ord = [0.1 0.5 1];
NS  = [1 5]./40;
Nc  = [65];
dBm = [1];

for k = 1:length(ord)
    for j = 1:length(NS)
        for i = 1:length(Nc)

            Coeff_vs_filt(NS(j),Nc(i),dBm(1) , symbols(1), n_prop_steps,etasp,ord(k));

        end
    end
   
end
