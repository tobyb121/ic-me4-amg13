clear all;
close all;

amg_cycle('v1',3,'v2',10,'v3',3,'smoother',@Jacobi);
k=500;
rv=zeros(k,7);
for i=1:7
    
    load('PoissonTestMatrices/poisson_80');
    A=-A;
    b=-b;
    x=ones(N,1);
       
    amg_cycle('reset');
    
    xv=x;
    WU=0;
    fprintf('Iterating:  setup');
    for j=1:k
        [xv,WUv]=amg_cycle(A,b,xv,1,i);
        WU=WUv+WU;
        rv(j,i)=norm(b-A*xv);
        fprintf('\b\b\b\b\b\b% 5d\n',j);
    end
end

semilogy(rv);