clear;

grid_rows=63;
grid_cols=63;

N=grid_rows*grid_cols;

B=ones(N,1)*[-1,-1,4,-1,-1];
B(grid_cols+1:grid_cols:N,4)=0;
B(grid_cols:grid_cols:N,2)=0;
d=[-grid_cols,-1,0,1,grid_cols];

A=spdiags(B,d,N,N);

clear('B','d');
X=peaks(grid_rows);
u=reshape(X,N,1);
b=A*u;
x=zeros(N,1);

xj=Jacobi(A,b,x,2/3,20);
x2=twoGridCorrection(A,x,b,grid_rows,grid_cols);
xv=vcycle(A,x,b,grid_rows,grid_cols,4);

Xj=reshape(xj,grid_rows,grid_cols);
X2=reshape(x2,grid_rows,grid_cols);
Xv=reshape(xv,grid_rows,grid_cols);

subplot(2,4,1);
surf(X);
subplot(2,4,2);
surf(Xj);
subplot(2,4,3);
surf(X2);
subplot(2,4,4);
surf(Xv);

subplot(2,4,6);
surf(Xj-X);
subplot(2,4,7);
surf(X2-X);
subplot(2,4,8);
surf(Xv-X);
