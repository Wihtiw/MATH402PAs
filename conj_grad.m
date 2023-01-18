% Author: Arda Can Türkan
%
% Description:
% Iterates using CG 
%
% A: matrix
% b: vector
% x0: initial vector
% maxit: maximum number of iterations
% tol: toleranca value
%
% Usage:
% [X, res, ite] = conj_grad(A,b,x0,tol,maxit)


function [X,res,ite] = conj_grad(A,b,x0,tol,maxit)
    x_k = x0;
    k=0;
    X = [x_k];
    r_k = A*x_k - b;
    p_k = -r_k;
    res = [norm(r_k)];
    while ((norm(r_k)>=tol) & k<=maxit)
        alpha_k = transpose(r_k)*r_k/(transpose(p_k)*A*p_k);
        x_k = x_k + alpha_k*p_k;
        X = [X,x_k];
        r_old = r_k;
        r_k = r_k + alpha_k*A*p_k;
        res = [res, norm(r_k)];
        beta_k = transpose(r_k)*r_k/(transpose(r_old)*r_old);
        p_k = -r_k + beta_k*p_k;
        k = k+1;
    end
    ite = k;
end
   