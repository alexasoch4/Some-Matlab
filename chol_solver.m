function [x] = chol_solver(A,b)
% Usage: [x] = chol_solver(A,b)
% Solve Ax = b by cholesky factorization
% Input:
% A = positive definite matrix
% b = right hand side
% Output:
% x = solution vector

L = chol(A,'lower'); % Matlab command
x = backsub(L', forsub(L, b));





