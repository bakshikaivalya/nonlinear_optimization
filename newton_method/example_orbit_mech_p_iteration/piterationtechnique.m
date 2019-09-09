function [p, v1_X, v2_X] = piterationtechnique(r1_X, r2_X, t1, t2, orbitdirection, tol, kmax)
% Author: Kaivalya Bakshi
% Note - This program assumed the use of earth canonical units for inputs,
% outputs and formulae used
% orbitdirection = 1 => short way orbit & orbidirection = 2 => long way
% orbit

% Computation of change in true anamoly depending on whether the orbit has
% the orbiting object moving in the 'short way' or the 'long way'

mu = 1;

% Given data
r1 = norm(r1_X);
r2 = norm(r2_X);

deltat = t2 - t1;

% Computation of deltaf depending on whether the orbit is short way or long
% way
    % Short way orbit
if (orbitdirection == 1)
    
    deltaf = acos(dot(r1_X, r2_X)/(norm(r1_X)*norm(r2_X)));
    
end

    % Long way orbit
if (orbitdirection == 2)
    
    deltaf = 2*pi - acos(dot(r1_X, r2_X)/(norm(r1_X)*norm(r2_X)));
    
end

% Definition of constants used in the program
k = r1*r2*(1 - cos(deltaf));
l = r1 + r2;
m = r1*r2*(1 + cos(deltaf));

% Calculation of p1, p2
p1 = k/(l + sqrt(2*m));

p2 = k/(l - sqrt(2*m));

% Choosing initial guess for p, viz. p0
p0 = (p1 + p2)/2;
    % For first iteration
p = p0;


% First iteration of Newton's Method
    % Computation of semi major axis for the first iteration
a = k/((2*m - l^2)*p^2 + 2*k*l*p - k^2);

    % Computation of F, G, Fdot
F = 1 - r2/p*(1 - cos(deltaf));

G = r1*r2*sin(deltaf)/sqrt(mu*p);

Fdot = sqrt(mu/p)*tan(deltaf)*((1 - cos(deltaf))/p - 1/r1 - 1/r2);

    % Computation of eccentric anamoly or hyperbolic anamloy depending on
    % whether the orbit is elliptic or hyperbolic for first iteration
if (a > 0)
    
   deltaE = atan2(-r1*r2*Fdot/sqrt(mu*a), 1 - r1/a*(1 - F));
   
   func0 = G + sqrt(a^3/mu)*(deltaE - sin(deltaE)) - deltat;
   dg_dp0 = -G/2/p - 3/2*a*(deltat - G)*(k^2 + (2*m - l^2)*p^2)/(m*k*p^2) + sqrt(a^3/mu)*(2*k*sin(deltaE))/(p*(k - l*p));
   
end

if (a < 0)
    
    deltaF = acosh(1 - r1/a*(1 - F));
    
    func0 = G + sqrt(-a^3/mu)*(sinh(deltaF) - deltaF) - deltat;
    dg_dp0 = -G/2/p - 3/2*a*(deltat - G)*(k^2 + (2*m - l^2)*p^2)/(m*k*p^2) + sqrt((-a)^3/mu)*(2*k*sin(deltaF))/(p*(k - l*p));
    
end

    % First iteration using Newton's Method
p(1) = p0 - func0/dg_dp0;
error(1) = abs(p(1) - p0);
pvalue = p(1);

% Applying Newton's Method for iterations beyond the first iteration
ind = 2;
while (error(ind-1) >= tol) && (ind <= kmax)
    
    % Computation of semi major axis for the kth iteration
    a = m*k*p/((2*m - l^2)*pvalue^2 + 2*k*l*pvalue - k^2);

    % Computation of F, G, Fdot
    F = 1 - r2/pvalue*(1 - cos(deltaf));

    G = r1*r2*sin(deltaf)/sqrt(mu*pvalue);

    Fdot = sqrt(mu/pvalue)*tan(deltaf)*((1 - cos(deltaf))/pvalue - 1/r1 - 1/r2);

    % Computation of eccentric anamoly or hyperbolic anamloy depending on
    % whether the orbit is elliptic or hyperbolic
    if (a > 0)
        
        deltaE = atan2(-r1*r2*Fdot/sqrt(mu*a), 1 - r1/a*(1 - F));
   
        func = G + sqrt(a^3/mu)*(deltaE - sin(deltaE)) - deltat;
        dg_dp = -G/2/pvalue - 3*a/2*(deltat - G)*(k^2 + (2*m - l^2)*pvalue^2)/(m*k*pvalue^2) + sqrt(a^3/mu)*(2*k*sin(deltaE))/(pvalue*(k - l*pvalue));
   
    end

    if (a < 0)
    
        deltaF = acosh(1 - r1/a*(1 - F));
    
        func = G + sqrt(-a^3/mu)*(sinh(deltaF) - deltaF) - deltat;
        dg_dp = -G/2/pvalue - 3*a/2*(deltat - G)*(k^2 + (2*m - l^2)*pvalue^2)/(m*k*pvalue^2) + sqrt((-a)^3/mu)*(2*k*sin(deltaF))/(pvalue*(k - l*pvalue));
    
    end
    
    % Calculation of value of g at x = x_(k - 1)
    g(ind - 1) = func;
    % Calculation of partial Derivative of g w.r.t x at x = x_(k - 1)
    dg_dp(ind - 1) = dg_dp;
    p(ind) = p(ind-1) - g(ind - 1)/dg_dp(ind - 1)
    error(ind) = abs(p(ind) - p(ind - 1));
    
    pvalue = p(ind);
    
    % Incrementing iteration index
    ind = ind+1;
            
end

Root = p(ind - 1);
p = Root;

% Computation of Gdot
Gdot = 1 - r1/p*(1 - cos(deltaf));

% Checking conservation of angular momentum holds
Residual = (F*Gdot - Fdot*G) - 1;

% Computation of v1 and v2
v1_X = (r2_X - F*r1_X)/G;

v2_X = Fdot*r1_X + Gdot*v1_X;

end