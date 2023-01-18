% Author: Arda Can Türkan
%
% Description:
% Iterates using Newton's Method with
% Armijo condition, stops if 
%
% Input:
% fhandle: function you like to use
% x0: initial x
% tol: tolerance
% maxit: max iterations
% alpha0: initial step length
% c: armijo constant
% mu: backtracking parameter
% amax: max number of armijo iteration
%
% Output:
% X: Matrix of past values of x
% Grad: Matrix of past values f x
% ite: number of iterations
%
% Usage:
% fhandle = @Rosenbrock
% x0 = [-0.5 1]
% tol = 1e-6
% maxit = 10000
% alpha0 = 1
% c = 1e-4
% mu = 0.5
% amax = 100
% [X, Grad, ite] = Newton_armijo(fhandle,x0,tol,maxit,alpha0,c,mu,amax)
%

function [X, Grad, ite] = Newton_armijo(fhandle,x0,tol,maxit,alpha0,c,mu,amax)
    X = [x0];
    ite = 0;
%    grad = gradient(fhandle);
%    hess = hessian(fhandle);
    [f, df, Hess] = fhandle(x0);

    k = 0;
    x_k = x0;
    Grad = [norm(df)];

    while ((norm(df) >= tol) & (k<= maxit))
        [f, df, Hess] = fhandle(x_k);
        pk = -(df / Hess);
        [alpha] = armijoalpha(fhandle, x_k, pk, alpha0, c, mu, amax);

        x_k = x_k + alpha*pk;
        [f, df , Hess] = fhandle(x_k);
        k = k+1;
        X = [X; x_k];
        Grad = [Grad; [norm(df)]];
    end
    ite = k;
end



