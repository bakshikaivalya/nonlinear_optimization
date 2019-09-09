%% Augmented Lagrangian Method implementation for unconstrained minimization
% of f: R^n -> R subject to the single inequality constraint g : R^n -> R
% Author: Kaivalya Bakshi

function [x_opt] = funcAugmentedLagrangeOpt(func, ineq_constr, x0, r_0_dash, ...
    lambda_0, gamma, max_iter)
% f = @func;
% g = @ineq_constr;

% Initializing the augmented Lagrangian methjod
p = 1;
r_p_dash = r_0_dash;
lambda_p = lambda_0;
x = x0;

% Executing the method for r_p incremented over maximum iterations
for(p = 1:1:max_iter)
   
    Psi = @(var) max(ineq_constr(var),-lambda_p/2/r_p_dash);
    
    A = @(var) func(var) + lambda_p*Psi(var) + r_p_dash*Psi(var)^2;
    
    [x_opt, x_opt_traj] = funcBFGSopt(A, x, max_iter);
    
    x = x_opt;
    
    lambda_p = lambda_p + 2*r_p*Psi(x);
    r_p_dash = gamma*r_p;
    
    p = p+1;
    
end


end