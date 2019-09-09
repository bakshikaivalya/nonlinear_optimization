function [rt_X, vt_X] = UniversalSingleToFSolverusingDiffipNewtonMethod(r0_X, v0_X, deltat, x0, tol, kmax)
% Function giving the position and velocity vectors of object in orbit for
% two body problem at time t = t0 + deltat after time of flight, deltat, 
% from starting time t0, given deltat and the position and velocity vectors
%  viz. r0_X and v0_X at t0 expressed in the same coordinate system X.
% Author: Kaivalya Bakshi
% Function Aguments-
% r0_X = Position vector of orbiting object w.r.t. to focus expressed in
% any orthonormal coordinate system X
% v0_X = Velocity of orbiting obect calculated using perifocal/ECI/Inertial
% reference frame expressed in any orthonormal coordinate system X 
% deltat = time of flight (s)
% x0 = Initial guess of root of variable x provided by the user. Do not use
% x0 = 0 as that will give an error since x0 appears in the denominator
% in expressions in the computations below and this would cause those 
% expressions to tend to infinity
% tol = Convergence criterion based on maximum permissible difference
% between values of root obtained in successive iterations
% kmax = Convergence criterion based on maximum number of permitted
% iterations

% Note- All the physical quantities in this program are in Earth Canonical
% Units

% Note- X must be the same for both r0_X and v0_X

% Note- This UniversalSingleToFSolver uses the r(x) expression to obtain
% the derivative of g(x) or g(x) - deltat the root for which we are trying
% to obtain

% Note- Universal time of flight problem algorithm cannot be used for
% parabolic orbits unless the infinite series formulae for C(z) and S(z)
% are used which is computatinally difficul and has not been executed 
% in this program. This program checks and prints if that is the case,
% so that in that case, program output is not useful

r0_xdotv0_X = dot(r0_X, v0_X);
r0 = norm(r0_X);

% Noting that we are using Earth Canonical Units
mu = 1;

% Obtaining value of a
h_X = cross(r0_X,v0_X);
e_X = 1/mu*(cross(v0_X,h_X) - mu*r0_X/r0);
h = norm(h_X);
e = norm(e_X);
p = h^2/mu;
a = p/(1 - e^2); % Note that a can also be calculated using the vis viva
                    % energy equation
                    
% Informing user if orbit is parabolic, then program output is not useful
% and this particular program for Universal time of flight algorithm cannot
% be used
    % Parabolic orbit =>
if (e == 1)
   fprintf('Orbit is parabolic and Universal Time of Flight algorithm cannot be used. Therefore please ignore program output since is not useful.');
end

str1 = num2str(r0_xdotv0_X);
str2 = num2str(r0);
str3 = num2str(mu/mu);
str4 = num2str(deltat);
str5 = num2str(a);

% Definition of functions C(z) and S(z) used as strings for elliptic and
% hyperbolic orbits
    % Elliptic orbit =>
if (e < 1)
    fprintf('Orbit is elliptic.');
    C = strcat('(1 - cos(sqrt(x^2/', str5, ')))/(x^2/', str5, ')');
    S = strcat('(sqrt(x^2/', str5, ')- sin(sqrt(x^2/', str5, ')))/', 'sqrt((x^2/', str5, ')^3)');
end
    % Hyperblic orbit =>
if (1 < e)
    fprintf('Orbit is hyperbolic.');
    C = strcat('(1 - cosh(sqrt(-x^2/', str5, ')))/(x^2/', str5, ')');
    S = strcat('-(sqrt(-x^2/', str5, ')- sinh(sqrt(-x^2/', str5, ')))/', 'sqrt((-x^2/', str5, ')^3)');
end

% Expression of function g(x) - deltat and dg(x)/dt as a string used as 
% argument to use in Newton's Method to calculate x satisfying g(x) = t(x) 
% - deltat = 0, i.e. the required root of g(x) assuming that mu = 1
str_gminusdeltat = strcat('(x^3*(', S, ') + ', str1, '/', str3, '*x^2*(', C, ') + ',...
    str2, '*x*', '(1 - x^2/', str5, '*(', S, ')))/', str3, '-', str4);

str_dg_dt = strcat('(x^2*(', C, ') + ', str1, '/', str3, '*x*(1 - x^2/', str5, '*(', S, ')) + ',...
    str2, '*(1 - x^2/', str5, '*(', C, ')))/', str3);

% Applying Newton's Method which requires inputting expression for the
% function g(x) - deltat the root of which we are finding as well as its
% derivative w.r.t the independent variable x 
func = str_gminusdeltat;
dfunc = str_dg_dt;

x = NewtonMethodFuncandDiffip(func, dfunc, x0, tol, kmax)

% Evaluation of functions C(z) and S(z) at given x for elliptic and
% hyperbolic orbits
    % Elliptic orbit =>
if (e < 1)
    C = (1 - cos(sqrt(x^2/a)))/(x^2/a);
    S = (sqrt(x^2/a) - sin(sqrt(x^2/a)))/sqrt((x^2/a)^3);
end
    % Hyperbolic orbit =>
if (1 < e)
    C = (1 - cosh(-sqrt(-x^2/a)))/(x^2/a);
    S = -(sqrt(-x^2/a) - sinh(sqrt(-x^2/a)))/sqrt((-x^2/a)^3);
end

% Magnitude of position vector at time t = t0 + deltat
r = x^2*C + r0_xdotv0_X*x*(1 - x^2/a*S)/sqrt(mu) + r0*(1 - x^2/a*C);

    % Calculating value of F, G, Fdot, Gdot
F = 1 - a/r0*(1 - cos(x/sqrt(a)));
G = deltat - x^3/sqrt(mu)*S;
Fdot = -sqrt(mu*a)/r/r0*sin(x/sqrt(a));
Gdot = 1 - x^2/r*C;

F_G_Fdot_Gdot = [F G Fdot Gdot];

rt_X = F*r0_X + G*v0_X;

vt_X = Fdot*r0_X + Gdot*v0_X;
   
end