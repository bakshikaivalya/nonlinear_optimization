function [c,ceq] = funcConstraint_eqns(x)
% AE 8803 (GER) Project
% Function returning constraint equations to use for nonlcon input argument
% to MATLAB's ga() function.


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


% Design state vector: x = [h b R r]

h = x(1);
b = x(2);
R = x(3);
r = x(4);

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

c(4) = r - R;

c(5) = omega - omega_max;

c(6) = -T + 32000; 


ceq = [];


end