function [x_opt, f_opt, fmin, fmax, favg, avgdist] = KaivalyaBakshi_GA(crossoverfraction, populationsize, generations, f)
% Genetic algorithm for AE 8803 GER deign project which gives output of
% optimal design state, objective function at optimal design state,
% objective function value bounds of min, max and the average for all the
% generations


global omega_max
global E
global S_e
global rho_s
global FOS

% rng(12345);

% Specifying design constants
omega_max = 10*pi; % -> rad/s % Maximum rotational speed of drive wheel
E = 29e6; % ->psi % Young's/Elastic modulus
S_e = 28e3*(32.18*12); % ->psi*(32.18*12)->lb/(in*s^2) % Endurance limit
rho_s = 0.28; % ->lb/in^3
FOS = 2; % Factor of safety


% Objective function to be minimized which includes application of the
% constraint equations implicitly
% f = @ funcObjective;


% GA parameters used
    % No elitism used
% crossoverfraction = 0.75;

% populationsize = 50;

% generations = 50;

mutationrate = 0.01;

% Initialization
    % We use the lower bounds and upper bounds specified to create random
    % individuals which lie within these bounds by using the MATLAB rand
    % function to generate a Gaussian distributed random variable which
    % lies in the closed interval [0,1]
    
    % Lower and upper bounds for the design variables
LB = [1 1 20 1]; UB = [10 10 60 60];

x(1:populationsize,1)=  ones(populationsize,1)*LB(1) ...
    + rand(populationsize,1)*(UB(1) - LB(1));
x(1:populationsize,2)=  ones(populationsize,1)*LB(2) ...
    + rand(populationsize,1)*(UB(2) - LB(2));
x(1:populationsize,3)=  ones(populationsize,1)*LB(3) ...
    + rand(populationsize,1)*(UB(3) - LB(3));
x(1:populationsize,4)=  ones(populationsize,1)*LB(4) ...
    + rand(populationsize,1)*(UB(4) - LB(4));

    % Objective function computation and ranking of the initial population
    % created
obj_fun=feval(f,x);
[obj_fun,j]=sort(obj_fun); x=x(j,:);


% Running the selection, cross over and mutation in loop to iterate
% through all the generations and arrive at the best design
    % Computation of variables needed inside the GA loop
populationfraction = round(crossoverfraction*populationsize);
w = flipud([1:populationfraction]'/sum([1:populationfraction]));
number_of_offspring = round((populationsize-populationfraction)/2);
for i = 1:1:generations
    
    
% Selection by Roulette selection
    % probability distribution function used for roulette selection
PDF = [0 cumsum(w(1:populationfraction))'];

    % Generating two random numbers
randomnumber1 = rand(1,number_of_offspring);
randomnumber2 = rand(1,number_of_offspring);
    
for k = 1:1:number_of_offspring
    for l = 2:1:populationfraction+1
        if randomnumber1(k)<=PDF(l) & randomnumber1(k)>PDF(l-1)
        A_index(k)=l-1;
        end
        if randomnumber2(k)<=PDF(l) & randomnumber2(k)>PDF(l-1)
        B_index(k)=l-1;
        end
    end       
k = k+1;    
end

% Crossover by Single Point Crossover
A_indices=1:2:populationfraction;
crossover_place=ceil(rand(1,number_of_offspring)*4);

for k=1:number_of_offspring
parents_concatenate=x(A_index(k),crossover_place(k))-x(B_index(k),crossover_place(k));

x(populationfraction+A_indices(k),:) = x(A_index(k),:);
x(populationfraction+A_indices(k)+1,:) = x(B_index(k),:);

x(populationfraction+A_indices(k),crossover_place(k))=x(A_index(k),crossover_place(k))- parents_concatenate;
x(populationfraction+A_indices(k)+1,crossover_place(k))=x(B_index(k),crossover_place(k))+ parents_concatenate;

if crossover_place(k)<4 % crossover when last variable not selected
x(populationfraction+A_indices(k),:) = [x(populationfraction+A_indices(k),1:crossover_place(k)) ...
    x(populationfraction+A_indices(k)+1,crossover_place(k)+1:4)];
x(populationfraction+A_indices(k)+1,:)=[x(populationfraction+A_indices(k)+1,1:crossover_place(k)) ...
    x(populationfraction+A_indices(k),crossover_place(k)+1:4)];
end % if
end

% Mutation by Uniform Mutation
MutatedRow=sort(ceil(rand(1,round((populationsize-1)*4*mutationrate))*(populationsize-1))+1);
MutatedCol=ceil(rand(1,round((populationsize-1)*4*mutationrate))*4);

for l=1:round((populationsize-1)*4*mutationrate)
LB_mutation = 1; UB_mutation = 10;
x(MutatedRow(l),MutatedCol(l))=(UB_mutation-LB_mutation)*rand + LB_mutation;
end

% Computation of objective function for mutated offspring chromosomes
obj_fun = feval(f,x);
[obj_fun,j]=sort(obj_fun); x=x(j,:);

favg(i+1) = mean(obj_fun);
fmin(i+1) = -(min(obj_fun) - mean(obj_fun));
fmax(i+1) = (max(obj_fun)- mean(obj_fun));

i = i+1;


end

% errorbar(1:1:generations+1,favg,fmin,fmax,'-');
% xlabel('Generations');
% ylabel('f');

x_opt = x(1,:);

obj_fun = @(x) (2*(2*x(3)+12)*x(1)*x(2) + 4*0.4*4*pi*x(3)^2);
f_opt = obj_fun(x_opt);

avgdist = pdist(x);


end