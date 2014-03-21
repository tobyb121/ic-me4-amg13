clear all;disp('Clearing');
close all;

disp('Setup');

k=200;

rows=64;
N=rows^2;

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

u=rand(N,1);

b=A*u;
x=zeros(N,1);

xk0=x;
rk0=norm(A*xk0-b);
%%
amg_cycle('v1',3,'v2',10,'v3',3,'smoother',@Jacobi);
xv=xk0;
rv=[rk0];
WU=0;
disp('Starting AMG cycle');
for i=1:k
    [xv,WUv]=amg_cycle(A,b,xv,1,9);
    WU=WUv+WU;
    rv=[rv,norm(A*xv-b)];
    disp(i);
end

[xj,rj]=Jacobi(A,b,x,round(WU));

xkn=xv;
rkn=norm(A*xkn-b);
rjn=norm(A*xj-b);

semilogy(linspace(0,length(rj),length(rj)),rj,linspace(0,length(rj),length(rv)),rv);

