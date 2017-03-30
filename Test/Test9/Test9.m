symbols      = [2^18];
n_prop_steps = 50;

Fn       = [7];
etasp    = [0.5 .*10.^(Fn/10)];
gamma = 1.2e-3;

NS  = [5,5]./40;
Nc  = [9,17,33,65];
dBm = [1,6];

for j = 1:length(symbols)
    for j = 1:length(NS)
        for i = 1:length(Nc)

            Coeff_vs_filt(NS(j),Nc(i),dBm(j) , symbols(1), n_prop_steps,etasp);

        end
    end
   
end
