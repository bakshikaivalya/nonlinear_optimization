% BFGS unconstrained optimization function to minimize a scalar function using the GSM method for line search
% Author: Kaivalya Bakshi
function [x_opt, x_opt_traj] = funcBFGSoptimization(func, x0, max_iter)

% Inputs- func = scalar function F : R^n -> R to be minimized, initial input
% value of x, max_iter = maximum number of iterations allowed

% Outputs- x_opt = optimal x which minimizes the function
                
[N,X] = size(x0);
% initializing BFGS Hessian approximation
H_tilde = eye(N,N);

% F = @func;

% bGoOn = true;
% Iters = 0;

% Computing initial search direction
eps = 0.001;
gradF_1 =  func_Gradient(func,x0,eps)
% grad1 = grad1';

iter = 0;
x = x0;

while(iter < max_iter)

    iter = iter + 1;
    
    % Obtaining search direction for present iteration k
    % S = S_k = -H_tilde_k * gradF_k-1
    S = -H_tilde * gradF_1;
    % normalize vector S
    S = S/norm(S);

    % Doing line search in the obtained search direction
    %   lambda = 1;
    %   lambda = linsearch(X, N, lambda, S, myFx);
    [x_opt_LineSearch, alpha_opt] = funcGoldenSearchMethod(func,x,gradF_1,1e-4);
    
    % Computing x_k using computed optimal alpha obtained from GSM
    % linesearch. This is the actual design state optimization update
    
    % d = lambda * S;
    % finding optimization step for present iteration, p = p_k = x_k - x_k-1
    p = alpha_opt * S;
    % finding x_k and assigning it to be the optimal design state at the
    % present iteration before storing the earlier optimizing design state
    % in the optimizer trajectory
    x_opt_traj(:,iter) = x;
    x = x + p;
    x_opt = x;
    
    % get gradient at x_k
    gradF_2 =  func_Gradient(func,x,eps)

    %   grad2 = grad2';
    % finding y = y_k = gradF_k - gradF_k-1
    y = gradF_2 - gradF_1
    % reassigning variables
    gradF_1 = gradF_2;

    % Tests for convergence
%     if (func(x_opt_traj(:,iter) - func(x)) < 0.0001)
%         return
%     end
    
    
%     for i = 1:N
%         if abs(d(i)) > DxToler(i)
%         break
%         end
%     end

%     if norm(gradF1) < gradToler
%         break
%     end

    % Computing BFGS approzimated Hessian for next iteration
    sigma = p'*y
    tau = y'*H_tilde*y;
    D = ((sigma + tau)/sigma^2)*(p*p') - (1/sigma)*(H_tilde*y*p' + p*(H_tilde*y)');
    H_tilde = H_tilde + D;
  
    
end


end                                     