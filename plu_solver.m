function [x]=plu_solver(A,b)
% Usage: [x]=lu_solver(A,b)
% Solves the linear system Ax = b by LU factorization LUx = b
% Input:
% A = matrix
% b = right hand side vector
% Output:
% x = solution vector
[L,U,P]=lu_fac_pivot(A) % Calculate LU factorization of a
x = backsub(U,forsub(L,P*b)); % Solve LUx = b