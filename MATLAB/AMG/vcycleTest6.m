clear all;
close all;

G={'uniform','graded'};

k=250;

zero=zeros(k,length(G));
ravj=zero;
rvj=zero;
rj=zero;
rgs=zero;
rcg=zero;

for i=1:length(G)
    
    load(['ADE_80_',G{i}]);
    %A=-A;
    %b=-b;
    x=ones(N,1);
       
    amg_cycle('reset');
    
    % AMG Jacobi
    amg_cycle('v1',3,'v2',10,'v3',3,'smoother',@Jacobi);
    xavj=x;
    WU=0;
    disp('AMG Jacobi');
    fprintf('Iterating:  setup');
    for j=1:k
        [xavj,WUv]=amg_cycle(A,b,xavj,1,10);
        WU=WUv+WU;
        ravj(j,i)=norm(b-A*xavj);
        fprintf('\b\b\b\b\b\b% 5d\n',j);
    end

    vcycle('v1',3,'v2',10,'v3',3,'smoother',@Jacobi);
    xvj=x;
    disp('Geometric V Jacobi');
    fprintf('Iterating:  setup');
    for j=1:k
        [xvj]=vcycle(A,b,xvj,80,80,6);
        rvj(j,i)=norm(b-A*xvj);
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
semilogy([ravj(:,1),rvj(:,1),rj(:,1),rcg(:,1)]/r0);
legend({'Algebraic Multigrid','Geometric Multgrid','Jacobi','Conjugate Gradient'});
figure;
semilogy([ravj(:,2),rvj(:,2),rj(:,2),rcg(:,2)]/r0);
legend({'Algebraic Multigrid','Geometric Multgrid','Jacobi','Conjugate Gradient'});
% figure(2);
% semilogy(rj);
% figure(3);
% semilogy(rcg);
% %%
% lrv=log10(rv);
% P=zeros(1,length(G));
% for i=1:length(G)
%    dn=[find(lrv(:,i)<-2,1,'first'),find(lrv(:,i)>-8,1,'last')];
%    P(i)=(dn(2)-dn(1))/10;
% end
% 
% figure(4);
% plot(G.^2,P);
% 
% %%
% lrj=log10(rj);
% P2=-(k-50)./(lrj(k,:)-lrj(50,:));
% lrcg=log10(rcg);
% P3=-(k-50)./(lrcg(k,:)-lrcg(50,:));
% 
% figure(5);
% plot(G.^2,P2);
% 
% figure(6);
% plot(G.^2,P3);
