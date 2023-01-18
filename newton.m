% Author: Arda Can Türkan
%
% Description:
% Iterates using Newton's Method until either 
% error is inside tolerance bounds or 
% maximum number of iterations has been reached.
%
% Input:
% f: Your function
% x_0: Initial guess
% tol: Tolerance bound, a positive real number. (eg. 1e-6)
% maxit: Maximum number of iterations.
%
% Output:
% x: Approximate solution to f(x) = 0
% hist: Value of x in each iteration
% hist_err: Value of error in each iteration
% iter: Number the function has iterated to get final x
%
%Usage:
% f(x) = (x^2+1)(x-1)
% x_0 = 1.8 initial guess
% tol = 1e-6
% maxit = 1000 number of max iterations
% [a, hist, hist_err, iter] = newton(f, 0, 1e-6, 1000)
%

function [x, hist, hist_err, iter] = Newton(f, x_0, tol, maxit)
    x = x_0;
    iter = 0;
    df = gradient(f);

    error = abs(f(x));
    hist_err(iter+1) = error;
    hist(iter+1) = x_0;

    while(abs(f(x)) >= tol)
        iter = iter + 1;

        der = df(x);
        x = x - f(x)/der;

        error = abs(f(x));
        hist(iter+1) = x;

        if(iter >= maxit)
            break;
        end
    end
    hist_err = abs(hist - x);
end
