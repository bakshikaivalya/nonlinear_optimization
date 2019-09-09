function [alpha_L, alpha_U] = func_Bracketing(func, x, s)

% Description: Bracketing for a nD function s.t. it is unimodal and 
% derivative of the function with search parameter alpha at initial point x
% is less than 0

% Inputs: @func = nD function taking vector x as input, x = nD initial
% point/vector for line search, s = search direction
% Outputs: Bracket or interval of alpha, [alpha_L, alpha_U] which contains
% the minimum


% Initializing
alpha_L = 0; alpha_U = 1; alpha1 = alpha_L; alpha3 = alpha_U;
Delta_alpha = alpha3 - alpha1;
s_cap = s/norm(s);

for(k = 1:1:100)
    
    
   x_alpha1 = x;
   x_alpha2 = x + Delta_alpha/2*s_cap;
   x_alpha3 = x + Delta_alpha*s_cap;
   
   f1 = func(x_alpha1);
   f2 = func(x_alpha2);
   f3 = func(x_alpha3);
   
   
   if(f3<min(f1,f2))
       
       Delta_alpha = 2*Delta_alpha; alpha3 = Delta_alpha;
       
   elseif(f1<min(f2,f3))
       
       
       Delta_alpha = Delta_alpha/2; alpha3 = Delta_alpha;
       
   else
       
       Delta_alpha = Delta_alpha; alpha3 = Delta_alpha;
       alpha_L = 0; alpha_U = alpha3;
       break;
    
    
   end


end


end