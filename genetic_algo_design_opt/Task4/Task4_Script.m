samples = 10;

f = @ funcObjective;
crossoverfraction = 0.75;
populationsize = 50;
generations = 50;

Davg = zeros(1,generations+1);
Fmin = zeros(1,generations+1); Fmax = zeros(1,generations+1);...
    Favg = zeros(1,generations+1);

for (i = 1:1:samples)
[x_opt, f_opt, fmin, fmax, favg, avgdist] = KaivalyaBakshi_GA(crossoverfraction, populationsize, generations, f);

Davg = Davg + avgdist';
Favg = Favg + favg;
Fmin = Fmin + fmin;
Fmax = Fmax + fmax;

i = i+1;
end

avgd = Davg/samples;
avgf = Favg/samples;
maxf = (Fmax-Favg)/samples;
minf = -(Fmin-Favg)/samples;

save ('KBakshi_GA.mat','avgd','avgf','maxf','minf')
