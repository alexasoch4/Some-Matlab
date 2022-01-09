%function needed for PLU_vs_Cholesky.m

function [x] = chol_solver(A,b)
L = chol(A,'lower'); 
x = backsub(L', forsub(L, b));





