% Author: Arda Can Türkan
%
% Description: Resonbrock banana function
%
% Input: 
% x: vector with 2 elements
%
% Output:
% f: value of the function
% df: Jacobian at the point x
% Hess: Hessian at the point x
%
% Usage:
% x = [-0.5 1]
% [f, df, Hess] = Rosenbrock(x)
%

function [f, df, Hess] = rosenbrock(x)
    x1 = x(1);
    x2 = x(2);
    f = 100*(x2-x1^2)^2 + (1-x1)^2;
    df = [400*x1^3 - 400*x1*x2 + 2*x1 - 2, -200*x1^2 + 200*x2];
    Hess = [1200*x1^2 - 400*x2 + 2, -400*x1; -400*x1, 200];
end