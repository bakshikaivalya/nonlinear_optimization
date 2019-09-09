% Continuous Genetic Algorithm
%
% minimizes the objective function designated in ff
% Before beginning, set all the parameters in parts
% I, II, and III
% Haupt & Haupt
% 2003
clc
clear all
close all
%_______________________________________________________
% I Setup the GA
ff=@ funcObjective; % objective function
npar=4; % number of optimization variables
% varhi=10; varlo=0; % variable limits
varhi(1)=10; varlo(1)=1; % variable limits
varhi(2)=10; varlo(2)=1; % variable limits
varhi(3)=60; varlo(3)=20; % variable limits
varhi(4)=60; varlo(4)=1; % variable limits
%_______________________________________________________
% II Stopping criteria
maxit=50; % max number of iterations
% mincost=-9999999; % minimum cost
%_______________________________________________________
% III GA parameters
popsize=50; % set population size
mutrate=0.01; % set mutation rate
crossoverfraction=0.75; % fraction of population kept
Nt=npar; % continuous parameter GA Nt=#variables
keep=round(crossoverfraction*popsize); % #population members that survive per generation update
nmut=round((popsize-1)*Nt*mutrate); % total number of mutations
M=round((popsize-keep)/2); % number of matings

% GA parameters used
    % No elitism used
crossoverfraction = 0.75;
mutationrate = 0.01;
populationsize = 50;

%_______________________________________________________
% Create the initial population
iga=0; % generation counter initialized
% par=(varhi-varlo)*rand(popsize,npar)+varlo; % random
% a =(varhi(1)-varlo(1))*rand(popsize,1)+varlo(1)*ones(popsize,1);
LB = [1 1 20 1]; UB = [10 10 60 60];

x(1:popsize,1)=  ones(populationsize,1)*LB(1) ...
    + rand(populationsize,1)*(UB(1) - LB(1));
x(1:popsize,2)=  ones(populationsize,1)*LB(2) ...
    + rand(populationsize,1)*(UB(2) - LB(2));
x(1:popsize,3)=  ones(populationsize,1)*LB(3) ...
    + rand(populationsize,1)*(UB(3) - LB(3));
x(1:popsize,4)=  ones(populationsize,1)*LB(4) ...
    + rand(populationsize,1)*(UB(4) - LB(4));
cost=feval(ff,x); % calculates population cost using ff
[cost,j]=sort(cost); % min cost in element 1
x=x(j,:); % sorted continuous least cost design state
% minc(1)=min(cost); % minc contains min of
% meanc(1)=mean(cost); % meanc contains mean of population
%_______________________________________________________
% Iterate through generations
M=round((popsize-keep)/2); % number of matings
prob=flipud([1:keep]'/sum([1:keep])); % weights chromosomes based on position in list
odds=[0 cumsum(prob(1:keep))']; % probability distribution function

for iga = 1:1:maxit
iga=iga+1; % increments generation counter
%_______________________________________________________
% Pair and mate

pick1=rand(1,M); % mate #1
pick2=rand(1,M); % mate #2
% ma and pa contain the indicies of the chromosomes
% that will mate
ic=1;
while ic<=M
for id=2:keep+1
if pick1(ic)<=odds(id) & pick1(ic)>odds(id-1)
ma(ic)=id-1;
end
if pick2(ic)<=odds(id) & pick2(ic)>odds(id-1)
pa(ic)=id-1;
end
end
ic=ic+1;
end
%_______________________________________________________
% Performs mating using single point crossover
ix=1:2:keep; % index of mate #1
xp=ceil(rand(1,M)*Nt); % crossover point
r=ones(1,M); % mixing parameter
for ic=1:M
xy=x(ma(ic),xp(ic))-x(pa(ic),xp(ic)); % ma and pa
% mate
x(keep+ix(ic),:)=x(ma(ic),:); % 1st offspring
x(keep+ix(ic)+1,:)=x(pa(ic),:); % 2nd offspring
x(keep+ix(ic),xp(ic))=x(ma(ic),xp(ic))-r(ic).*xy;
% 1st
x(keep+ix(ic)+1,xp(ic))=x(pa(ic),xp(ic))+r(ic).*xy;
% 2nd
if xp(ic)<npar % crossover when last variable not selected
x(keep+ix(ic),:)=[x(keep+ix(ic),1:xp(ic)) x(keep+ix(ic)+1,xp(ic)+1:npar)];
x(keep+ix(ic)+1,:)=[x(keep+ix(ic)+1,1:xp(ic)) x(keep+ix(ic),xp(ic)+1:npar)];
end % if
end
%_______________________________________________________
% Mutate the population
mrow=sort(ceil(rand(1,nmut)*(popsize-1))+1);
mcol=ceil(rand(1,nmut)*Nt);

for ii=1:nmut
% par(mrow(ii),mcol(ii))=(varhi-varlo)*rand+varlo;
varhiC = 10; varloC = 1;
x(mrow(ii),mcol(ii))=(varhiC-varloC)*rand+varloC;
% mutation
end % ii
%_______________________________________________________
% The new offspring and mutated chromosomes are evaluated
cost=feval(ff,x);
%_______________________________________________________
% Sort the costs and associated parameters
[cost,ind]=sort(cost);
x=x(ind,:);
%_______________________________________________________
% Do statistics for a single nonaveraging run
minc(iga+1)=min(cost);
meanc(iga+1)=mean(cost);
%_______________________________________________________
% Stopping criteria
if iga>maxit %| cost(1)<mincost
break
end
[iga cost(1)]
end %iga

% Displays the output
% day=clock;
% disp(datestr(datenum(day(1),day(2),day(3),day(4),day(5),...
% day(6)),0))
% disp(['optimized function is ' ff])
format short g
disp(['popsize = ' num2str(popsize) ' mutrate = '...
num2str(mutrate) ' # par = ' num2str(npar)])
disp(['#generations=' num2str(iga) ' best cost=' num2str(cost(1))])
disp(['best solution'])
disp([num2str(x(1,:))])
disp('continuous genetic algorithm')

figure(24)
iters=0:length(minc)-1;
plot(iters,minc,iters,meanc, '-');
xlabel('generation');ylabel('cost');
text(0,minc(1),'best'); text(1,minc(2),'population average')