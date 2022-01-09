%function needed for
function [x]=plu_solver(A,b)
[L,U,P]=lu_fac_pivot(A) 
x = backsub(U,forsub(L,P*b)); 
