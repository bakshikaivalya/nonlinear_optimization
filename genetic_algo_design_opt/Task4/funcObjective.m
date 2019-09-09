function [fitness] = funcObjective(state)
% Fittness = Volume = 28l*h*b + 4*0.4*4*pi*R^2

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

[m,n]=size(state);


for i=1:m
    

h=state(i,1);
% if h>10
% h=10;
% elseif h<1
% h=1;
% end    

b=state(i,2);
% if b>10
% b=10;
% elseif b<1
% b=1;
% end

R=state(i,3);
% if R>60
% R=60;
% elseif R<20
% R=20;
% end

r=state(i,4);
% if r>60
% r=60;
% elseif r<1
% r=1;
% end


% Design constraints
l = 2*R+12; % Length of connecting rod

% Area moments of inertia
I_h = b*h^3/12;
I_b = h*b^3/12;

T = 0.85*200*24^2*r/R;

f_h = pi^2*E*I_h/l^2;
c(1) = -f_h + FOS*R*T/(4*r);

f_b = pi^2*E*I_b/(2*l^2);
c(2) = -f_b + FOS*R*T/(4*r);

v_max = 70*17.6; % mph*17.6 -> in/s
omega = v_max/R;
a_c = omega^2*r; % Maximum centripetal acceleration of connecting rod
M = rho_s*l*h*b * l * a_c / 8; % Maximum bending moment of connecting rod
c(3) = FOS*M*h/(2*I_h) - S_e;

c(4) = omega - omega_max;

c(5) = -T + 32000;

c(6) = r - R;

c(7) = 1 - h;
c(8) = h - 10;

c(9) = 1 - b;
c(10) = b - 10;

c(11) = 20 - R;
c(12) = R - 60;

c(13) = 1 - r;
c(14) = r - 60;

fitness(i) = 2*(2*R+12)*h*b + 4*0.4*4*pi*R^2;

% Applying the constraint equations violation cost to the objective
% function evaluation
for j = 1:1:14
    if c(j)>0 && j<6
        fitness(i) = fitness(i) + 1e6*c(j)/abs(c(j));
    end
    if c(j)>0 && j>=6
        fitness(i) = fitness(i) + 1e6*c(j)/abs(c(j));    
    end
end


end


end