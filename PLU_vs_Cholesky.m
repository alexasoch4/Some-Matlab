%comparing PLU and Cholesky factorization for a given matrix A and b. Want to see which one is more efficient in time and has less error.

time_chol = [];
time_plu = [];
error_chol = [];
error_plu = [];

for k = 2:1:10 
    n = 2^k;
    h = 1/n; 
    v1 = ones(n-1,1); 
    v2 = ones(n-2,1); 
    
    %set up A
    A = diag(v1*(4*h^2+2)) + diag(v2*(-1), -1) + diag(v2*(-1),1);
    x = h:h:1-h;
    x = x';
    
    yexact = 1/4/(exp(2)-exp(-2))*(exp(2*x)-exp(-2*x))+1/4*x;
    %set up b
    b = x*h^2;
    b(1) = b(1) + 0;
    b(n-1) = b(n-1) + 1/2;
    b;
   
    %cholesky factorization
    tchol_start = tic;
    ychol = chol_solver(A,b);
    tchol_end = toc(tchol_start);
    time_chol = [time_chol, tchol_end];
    error_chol = [error_chol, norm(ychol-yexact,inf)];
    
    %plu factorization
    tplu_start = tic;
    yplu = plu_solver(A,b);
    tplu_end = toc(tplu_start);
    time_plu = [time_plu, tplu_end];
    error_plu = [error_plu, norm(yplu - yexact, inf)];
end 

k = 2:1:10;
h = 2.^(-k);

%elasped time comparison
figure 
plot(log(h), time_plu, 'r-*');
hold on 
plot(log(h), time_chol, 'b-o');
xlabel('log(h)');
ylabel('elapsed time');
legend('PLU', 'Cholesky');
hold off

%error comparison
figure 
plot(log(h), log(abs(error_plu)), 'r-*');
hold on 
plot(log(h), log(abs(error_chol)), 'b-o');
xlabel('log(h)');
ylabel('log(error)');
legend('PLU', 'Cholesky');
hold off


    
