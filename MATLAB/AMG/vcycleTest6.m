clear all;
close all;

G={'uniform','graded'};

k=250;

zero=zeros(k,length(G));
ramg=zero;
rvj=zero;
rj=zero;
rgs=zero;
rcg=zero;

for i=1:length(G)
    
    load(['ADE_80_',G{i}]);
    x=ones(N,1);
       
    amg_cycle('reset');
    
    % AMG
    amg_cycle('v1',3,'v2',10,'v3',3,'smoother',@Jacobi);
    xamg=x;
    WU=0;
    disp('AMG Jacobi');
    fprintf('Iterating:  setup');
    for j=1:k
        [xamg,WUv]=amg_cycle(A,b,xamg,1,10);
        WU=WUv+WU;
        ramg(j,i)=norm(b-A*xamg);
        fprintf('\b\b\b\b\b\b% 5d\n',j);
    end
    
    fprintf('Iterating: Jacobi %d\n',WU);
    [xj,rjn]=Jacobi(A,b,x,k*ceil(WU/k));
    rjn=rjn(1:length(rjn)/k:length(rjn));
    rj(:,i)=rjn;
    fprintf('Iterating: CG %d\n',WU);
    [xcg,rcgn]=conjgrad(A,b,x,k*ceil(WU/k));
    rcgn=rcgn(1:length(rcgn)/k:length(rcgn));
    rcg(:,i)=rcgn;
end
r0=norm(b-A*x);
figure;
semilogy([ramg(:,1),rj(:,1),rcg(:,1)]/r0);
legend({'Algebraic Multigrid','Jacobi','Conjugate Gradient'});
figure;
semilogy([ramg(:,2),rj(:,2),rcg(:,2)]/r0);
legend({'Algebraic Multigrid','Jacobi','Conjugate Gradient'});
