% Description: Golden section method to find the optimal x minimizing a
% function on R^n starting from a given initial point and a given search
% direction and the acceptable tolerance
% Author: Kaivalya Bakshi

function [x_opt_LineSearch, alpha_opt] = GSM_BakshiKaivalya(func,x,s,eps)

% Inputs:
%   func = function on R^n
%   x = initial point \in R^n to start line search from
%   eps = tolerance = ratio of final interval width to initial interval
%   width of the line search
% Outputs: 
%   x_Opt_LineSearch = Optimal x minimizing the function starting from x
%   along direction s


% Initializing
tau = 0.38197;
N = log(eps)/log(1 - tau) + 4;

[alpha_U, alpha_L] = funcBracketing(func,x,s);

s_cap = s/norm(s);

Delta_alpha = alpha_U - alpha_L;


% Initializing golden section method
tau = 0.38197;
alpha1 = (1 - tau)*alpha_L + tau*alpha_U;
alpha2 = tau*alpha_L + (1- tau)*alpha_U;

x_L = x; f_L = func(x_L);
x1 = x + alpha1*s_cap; f1 = func(x1);
x2 = x + alpha2*s_cap; f2 = func(x2);
x_U = x + Delta_alpha*s_cap; f_U = func(x_U);


% Storing alpha values for postprocessing
Alpha(1:4,1) = [alpha_L; alpha1; alpha2; alpha_U]; 


for(k = 1:1:floor(N))
    
   
    if(f1 < min(f_L, f2))
        
        % throw away current alpha_U
        % reassigned values
        % alpha_L = alpha_L;
        alpha_U = alpha2; x_U = x2; f_U = f2;
        % Delta_alpha = alpha_U - alpha_L;
        alpha2 = alpha1; x2 = x1; f2 = f1;
        
        % computed new point
        alpha1 = (1 - tau)*alpha_L + tau*alpha_U;
        x1 = x + alpha1*s_cap; f1 = func(x1);
        
        % Storing alpha values for postprocessing
        Alpha(4 + k,1) = alpha1;

        
    elseif(f2 < min(f1, f_U))
        
        % throw away current alpha_L
        % reassigned values
        % alpha_U = alpha_U;
        alpha_L = alpha1; f_L = f1;
        % Delta_alpha = alpha_U - alpha_L;
        alpha1 = alpha2; f1 = f2;
        
        % computed new point
        alpha2 = tau*alpha_L + (1- tau)*alpha_U;
        x2 = x + alpha2*s_cap; f2 = func(x2);
        
        % Storing alpha values for postprocessing
        Alpha(4 + k,1) = alpha2;
        
    end
    
    
end


x_opt_LineSearch = x + (alpha1 + alpha2)/2*s_cap;
alpha_opt = (alpha1 + alpha2)/2;
    
    
end