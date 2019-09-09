function [Root] = NewtonMethodFuncandDiffip(func, dfunc, x0, tol, kmax)
% Function that calculates roots for a real valued single variable function using Newton's
% Method
% Author: Kaivalya Bakshi
% Function arguments-
% func = expression for function of single variable in single quotes, ''
% dfunc = expression for dierivative of single variable w.r.t independent
% variable x in single quotes, ''
% x0 = Initial guess of root provided by the user
% tol = Convergence criterion based on maximum permissible difference
% between values of root obtained in successive iterations
% kmax = Convergence criterion based on maximum number of permitted
% iterations

% The function name is denoted by g in the algorithm as well as the data
% table

% Obtaining the expression of the real valued function and of the
% derivative of the real valued function w.r.t the independent variable
% from the function arguments 
func = inline(func);
dfunc = inline(dfunc);

% Applying Newton's Method for the first iterations
% Calculation of partial Derivative of g w.r.t x at x = x0
dg_dx0 = dfunc(x0);
x(1,1) = x0 - (func(x0)/dg_dx0);
error(1) = abs(x(1)-x0);

% Applying Newton's Method for iterations beyond the first iteration
k = 2;
while (error(k-1) >= tol) && (k <= kmax)
    
    % Calculation of value of g at x = x_(k - 1)
    g(1,k - 1) = func(x(1,k - 1));
    % Calculation of partial Derivative of g w.r.t x at x = x_(k - 1)
    dg_dx(1,k - 1) = dfunc(x(1,k - 1));
    x(1,k) = x(1,k-1) - func(x(1,k - 1))/dg_dx(1,k - 1);
    error(1,k) = abs(x(1,k) - x(1,k - 1));
    
    % Creating vectors for creating table of results
    g(1,k) = func(x(1,k));
    dg_dx(1,k) = dfunc(x(1,k));
    
    % Incrementing iteration index
    k = k+1;
        
end

Root = x(1,k - 1);

% Printing table of [Iteration number, x_k, g(x_k), g'(x_k), '|x_k - x_(k - 1)|']
Table = [[1:1:k - 1]', x', g', dg_dx', error'];
xlswrite('DataTable_AllIterations.xlsx', Table);

end