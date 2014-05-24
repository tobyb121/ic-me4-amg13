clear all;
close all;

amg_cycle('v1',3,'v2',10,'v3',3,'smoother',@Jacobi);

G={'graded','uniform'};

k=500;

rv=zeros(k,length(G));
rj=zeros(k,length(G));
rcg=zeros(k,length(G));
for i=1:length(G)
    
    load(['ADE_80_',G{i}]);
    %A=-A;
    %b=-b;
    x=ones(N,1);
       
    amg_cycle('reset');
    
    xv=x;
    WU=0;
    fprintf('Iterating:  setup');
    for j=1:k
        [xv,WUv]=amg_cycle(A,b,xv,1,10);
        WU=WUv+WU;
        rv(j,i)=norm(b-A*xv);
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

figure(1);
semilogy(rv);
figure(2);
semilogy(rj);
figure(3);
semilogy(rcg);
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
