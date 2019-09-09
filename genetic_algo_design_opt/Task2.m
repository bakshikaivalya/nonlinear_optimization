clc
clear all
% close all

rng(12345);

elitecount = 2;
crossoverfraction = 0.75;
mutationrate = 0.01;
populationsize = 50;

baseoptions = gaoptimset(@ga);
baseoptions = gaoptimset(baseoptions,'PlotFcns',{@gaplotdistance, @gaplotrange},...
'PlotInterval',1,...
'NonlinConAlgorithm','penalty',...
'EliteCount',elitecount,...
'SelectionFcn',@selectionroulette,...
'CrossoverFcn',@crossoverscattered,...
'CrossoverFraction',crossoverfraction,...
'MutationFcn',{@mutationuniform, mutationrate},...
'PopulationSize',populationsize,...
'Generations',50,...
'OutputFcns',@gaoutputfun);

global omega_max
global E
global S_e
global rho_s
global FOS


% Specifying design constants
omega_max = 10*pi; % -> rad/s % Maximum rotational speed of drive wheel
E = 29e6; % ->psi % Young's/Elastic modulus
S_e = 28e3*(32.18*12); % ->psi*(32.18*12)->lb/(in*s^2) % Endurance limit
rho_s = 0.28; % ->lb/in^3
FOS = 2; % Factor of safety


% Specifying the objective function
obj_fun = @(x) (2*(2*x(3)+12)*x(1)*x(2) + 4*0.4*4*pi*x(3)^2);
% obj_fun = @funcFitness;


%Specifying constraint equations
    % Specifying inequality constraints using another MATLAB function
nonlcon = @funcConstraint_eqns;
    % Side constraints
LB = [1, 1, 20, 1];
UB = [10, 10, 60, 60];

x = ga(obj_fun,4,[],[],[],[],LB,UB,nonlcon,baseoptions)

% i = 0; samples = 10;
% generations = 50;
% Davg = zeros(1,generations+1);
% Fmin = zeros(1,generations+1); Fmax = zeros(1,generations+1);...
%     Favg = zeros(1,generations+1);
% 
% for (i = 1:1:samples)
% x = ga(obj_fun,4,[],[],[],[],LB,UB,nonlcon,baseoptions);
% 
% Davg = Davg + avgdist;
% Favg = Favg + favg;
% Fmin = Fmin + fmin;
% Fmax = Fmax + fmax;
% 
% i = i+1;
% end
% 
% avgd = Davg/samples;
% avgf = Favg/samples;
% maxf = (Fmax-Favg)/samples;
% minf = -(Fmin-Favg)/samples;
% 
% save ('popsize100avg.mat','avgd','avgf','maxf','minf')
% 
% figure(1)
% subplot(2,1,1);
% xlim = [1,generations+1];
% plot(1:1:generations+1,avgd,'.');
% grid on
% title('Average distance averaged over 10 optimized samples')
% subplot(2,1,2);
% xlim = [1,generations+1];
% % plot(1:1:generations+1,avgf);
% % errorbar(1:1:generations+1,avgf,minf,maxf,'-');
% shadedErrorBar(1:1:generations+1, avgf, [maxf; minf],'-');
% grid on
% title('Average f averaged over 10 optimized samples')

