close all;
clear all;

grid_rows=20;
grid_cols=20;

N=grid_rows*grid_cols;

load('..\Matrices\GeometricTest_1_20131203\matrices.mat');

A=A6;
b=b6;

x=zeros(N,1);
%%
v1=10;
v2=10;
v3=40;
n=4;
vcycles=100;

WU=(v1+v3)*(4*(1-4^-n)/3)+v2*4^-n;
k=round(vcycles*WU);

disp(['Work Units: ',num2str(vcycles),'*',num2str(WU),'=',num2str(vcycles*WU)]);
%%

[xj,Rj]=Jacobi(A,b,x,k);

[xgs,Rgs]=GaussSeidel(A,b,x,k);

vcycle('v1',v1,'v2',v2,'v3',v3,'smoother',@Jacobi)
xvj=x;
Rvj=zeros(vcycles,1);
for i=1:vcycles
    xvj=vcycle(A,b,xvj,grid_rows,grid_cols,n);
    Rvj(i)=norm(A*xvj-b);
end

vcycle('v1',v1,'v2',v2,'v3',v3,'smoother',@GaussSeidel)
xvgs=x;
Rvgs=zeros(vcycles,1);
for i=1:vcycles
    xvgs=vcycle(A,b,xvgs,grid_rows,grid_cols,n);
    Rvgs(i)=norm(A*xvgs-b);
end

rj=A*xj-b;
rgs=A*xgs-b;
rv=A*xvj-b;
rv2=A*xvgs-b;

Xj=reshape(xj,grid_rows,grid_cols);
Xgs=reshape(xgs,grid_rows,grid_cols);
Xvj=reshape(xvj,grid_rows,grid_cols);
Xvgs=reshape(xvgs,grid_rows,grid_cols);

subplot(2,2,1);surf(Xj);axis tight
title(['Jacobi (',num2str(k),' iterations) ||r||=',num2str(norm(rj))]);
subplot(2,2,2);surf(Xvj);axis tight
title(['Multigrid V-Cycle (',num2str(n),' Grids, Jacobi) ||r||=',num2str(norm(rv))]);
subplot(2,2,3);surf(Xgs);axis tight
title(['Gauss Seidel (',num2str(k),' iterations) ||r||=',num2str(norm(rgs))]);
subplot(2,2,4);surf(Xvgs);axis tight
title(['Multigrid V-Cycle (',num2str(n),' Grids, Gauss Seidel) ||r||=',num2str(norm(rv2))]);

figure;
subplot(2,2,1);semilogy(Rj);axis tight
title(['Jacobi (',num2str(k),' iterations) ||r||=',num2str(norm(rj))]);
subplot(2,2,2);semilogy(Rvj);axis tight
title(['Multigrid V-Cycle (',num2str(n),' Grids, Jacobi) ||r||=',num2str(norm(rv))]);
subplot(2,2,3);semilogy(Rgs);axis tight
title(['Gauss Seidel (',num2str(k),' iterations) ||r||=',num2str(norm(rgs))]);
subplot(2,2,4);semilogy(Rvgs);axis tight
title(['Multigrid V-Cycle (',num2str(n),' Grids, Gauss Seidel) ||r||=',num2str(norm(rv2))]);
