clear;

k=10000;

grid_rows=20;
grid_cols=20;

N=grid_rows*grid_cols;

load('..\Matrices\GeometricTest_1_20131203\matrices.mat');

x=zeros(N,1);
%%
xj=Jacobi(A,b,x,k);
xgs=GaussSeidel(A,b,x,k);
vcycle('v1',10,'v2',10,'v3',20,'smoother',@Jacobi)
xv=vcycle(A,b,x,grid_rows,grid_cols,3);
vcycle('v1',10,'v2',10,'v3',20,'smoother',@GaussSeidel)
xv2=vcycle(A,b,x,grid_rows,grid_cols,3);

rj=A*xj-b;
rgs=A*xgs-b;
rv=A*xv-b;
rv2=A*xv2-b;

Xj=reshape(xj,grid_rows,grid_cols);
X2=reshape(xgs,grid_rows,grid_cols);
Xv=reshape(xv,grid_rows,grid_cols);
Xv2=reshape(xv2,grid_rows,grid_cols);

subplot(2,4,1);surf(Xj);axis tight
title(['Jacobi (',num2str(k),' iterations) ||r||=',num2str(norm(rj))]);
subplot(2,4,2);surf(X2);axis tight
title(['Gauss Seidel (',num2str(k),' iterations) ||r||=',num2str(norm(rgs))]);
subplot(2,4,3);surf(Xv);axis tight
title(['Multigrid V-Cycle (5 Grid) ||r||=',num2str(norm(rv))]);
subplot(2,4,4);surf(Xv2);axis tight
title(['Multigrid V-Cycle (5 Grid) ||r||=',num2str(norm(rv2))]);