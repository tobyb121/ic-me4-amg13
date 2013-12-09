clear all;
close all;

k=100;

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
%%
vcycleDbg('v1',3,'v2',10,'v3',3,'smoother',@GaussSeidel)
xv=vcycleDbg(A,b,x,grid_rows,grid_cols,3,x,b);

plot_vcycle;