% MATLAB script which plots deltat as a functoin of p in the p iteration
% technique assuming that the orbit is elliptic

% Note - The graphed domain must be manipulated manually after looking at
% some preliminary graphs so as to avoid
% having complex numbers passed as arguments to the atan2() function

mu = 1;

% Short way orbit
if (orbitdirection == 1)
    
    deltaf = acos(dot(r1_X, r2_X)/(norm(r1_X)*norm(r2_X)));
    
end

% Long way orbit
if (orbitdirection == 2)
    
    deltaf = 2*pi - acos(dot(r1_X, r2_X)/(norm(r1_X)*norm(r2_X)));
    
end

% Definition of constants used in the program
r1 = norm(r1_X);
r2 = norm(r2_X);

% Time interval
deltat = t2 - t1;

k = r1*r2*(1 - cos(deltaf));
l = r1 + r2;
m = r1*r2*(1 + cos(deltaf));

% Calculation of p1, p2
p1 = l/(l + sqrt(m))

p2 = 1/(1 - sqrt(m));

% Domain of graph
p = p1:0.01:40;

% Obtaining deltat as a function of p
    %Calculation of a as a function of p assuming the orbit is elliptic
a = k/((2*m^2 - l^2)).*p.^2 + 2*k*l.*p - k^2;

    % Computation of F, G, Fdot as functions of p
F = 1 - r2./p*(1 - cos(deltaf));

G = r1*r2*sin(deltaf)./sqrt(mu.*p);

Fdot = sqrt(mu./p)*tan(deltaf/2).*((1 - cos(deltaf)./p) - 1/r1 - 1/r2);

    % Calculation of deltaE as a function of p 
   deltaE = atan2(-r1*r2.*Fdot./sqrt(mu.*a), 1 - r1./a.*(1 - F));
   
   % Calculation of deltat as a function of p
deltat = G + sqrt(a.^3./mu).*(deltaE - sin(deltaE));

% Plotting deltat as a function of p
plot(p, deltat);
grid on;




