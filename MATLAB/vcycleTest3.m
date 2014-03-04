clear all;
close all;

k=10;

grid_rows=63;
grid_cols=63;

N=grid_rows*grid_cols;

B=ones(N,1)*[-1,-1,4,-1,-1];
B(grid_cols+1:grid_cols:N,4)=0;
B(grid_cols:grid_cols:N,2)=0;
d=[-grid_cols,-1,0,1,grid_cols];

A=spdiags(B,d,N,N);

clear('B','d');
%X=peaks(grid_rows);
X=rand(grid_rows,grid_cols);
X=X-mean(mean(X));
u=reshape(X,N,1);
b=A*u;
x=rand(N,1);
x=x-mean(x);

figure;
X0=reshape(x,grid_rows,grid_cols);
surf(abs(fft2(X0-X)));


xk0=x;
rk0=norm(A*xk0-b);
%%
vcycle('v1',2,'v2',5,'v3',2,'smoother',@Jacobi)
xv=vcycle(A,b,x,grid_rows,grid_cols,5);

Xv=reshape(xv,grid_rows,grid_cols);

xj=Jacobi(A,b,x,5);
Xj=reshape(xj,grid_rows,grid_cols);

figure;
surf(abs(fft2(Xv-X)));
zlim([0,50]);caxis([0,50]);
set(gca,'visible','off');

figure;
surf(abs(fft2(Xj-X)));
zlim([0,50]);caxis([0,50]);
set(gca,'visible','off');