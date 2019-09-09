function gradF = func_Gradient(func,x,eps)
%% Finds the gradient of a scalar function F: R^n -> R w.r.t its argument using forward finite difference method

% Inputs- func = function name of scaler function F saved as *.m file, 
% x \in R^n at which gradF is to be found, eps = finite difference step 
% size used

% Outputs- gradF = gradient of F at x

%Dimensions of function input x-
n = size(x);

F_of_x = func(x);

for(i = 1:1:n)
    
    x_s = x;
    x_s(i) = x(i) + eps;
    gradF(i,1) = (func(x_s) - F_of_x)/eps;    
    
end

gradF = gradF;


end