%from Q3

time_chol = [];
time_plu = [];
error_chol = [];
error_plu = [];

for k = 2:1:10 %change
    n = 2^k; %change
    h = 1/n; %change
    v1 = ones(n-1,1); %change
    v2 = ones(n-2,1); %change
    
    %set up A and b 
    %change
    A = diag(v1*(4*h^2+2)) + diag(v2*(-1), -1) + diag(v2*(-1),1);
    %change
    x = h:h:1-h;
    x = x';
    %change
    yexact = 1/4/(exp(2)-exp(-2))*(exp(2*x)-exp(-2*x))+1/4*x;
    %change
    b = x*h^2;
    b(1) = b(1) + 0;
    b(n-1) = b(n-1) + 1/2;
    b;
     

%     
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

%compare elasped time 
figure 
plot(log(h), time_plu, 'r-*');
hold on 
plot(log(h), time_chol, 'b-o');
xlabel('log(h)');
ylabel('elapsed time');
legend('PLU', 'Cholesky');
hold off

%compare error
figure 
plot(log(h), log(abs(error_plu)), 'r-*');
hold on 
plot(log(h), log(abs(error_chol)), 'b-o');
xlabel('log(h)');
ylabel('log(error)');
legend('PLU', 'Cholesky');
hold off


    