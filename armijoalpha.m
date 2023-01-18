% Author: Arda Can Türkan
%
% Description:
% a function that tries to reach a step-length which satisfies
% the Armijo condition for the given x
% 
% Input:
% fhandle : the function that we are operate on
% x0      : current iterate
% alpha0  : initial step length
% c       : Armijo constant
% mu      : backtracking parameter
% amax    : maximum number of iterations for Armijo condition
%
% Output:
% alpha: a suitable step length alpha 
%
% Usage:
% x = [-0.5 1]
% c = 1e-4
% mu = 0.5
% amax = 100
% alpha0 = 1
% call the armijo function as,
% [alpha] = armijo(@Rosenbrock, x, p, alpha0, c, mu, amax)

function [alpha] = armijoalpha(fhandle, x, p, alpha0, c, mu, amax)
    j = 0;
    alpha = alpha0;
    [f, df, Hess] = fhandle(x);
    while (fhandle(x + alpha * p) >= fhandle(x) + c * alpha * transpose(df) * p) & (amax >= j)
            alpha = alpha * mu;
            j = j + 1;
    end
end