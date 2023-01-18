% Author: Arda Can Türkan
%
% Description:
% Iterates using SR1 with
% Armijo condition
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
% H0 = I
% [X, Grad, ite] = SR1_inverse(fhandle,x0,tol,H0,maxit,alpha0,c,mu,amax)
%

function [X,Grad,ite]=SR1_inverse(fhandle,x0,tol,H0,maxit,alpha0,c,mu,amax)
    X = [x0];
    x_k = x0
    [f,df, Hess] = fhandle(x_k);
    Grad = [norm(df)];
    H_k = H0;
    k=0;
    while ((norm(df) >= tol) & (k<= maxit))
        [f,df,Hess] = fhandle(x_k);
        x_old = x_k;
        df_old = df;
        p_k = -H_k*df;
        [alpha] = armijoalpha(fhandle, x_k, p_k, alpha0, c, mu, amax);
        x_k = x_k + alpha*p_k;
        [f,df,Hess] = fhandle(x_k);
        X = [X,x_k];
        Grad = [Grad,norm(df)];
        s_k = x_k - x_old;
        y_k = df - df_old;
        H_k = H_k + (s_k - H_k*y_k)*transpose(s_k-H_k*y_k)/(transpose(s_k-H_k*y_k)*y_k);
        k = k+1;
    end
    ite = k;
    
end