%% vcycleTest
% Performs iterations of amg v-cycle on a 64x64 element grid using 2D
% laplace operator:
%       (0  -1  0)
%       (-1  4 -1)
%       (0  -1  0)
% Compares result to iteration of Jacobi (for same number of WU) and plots
% residuals for both cases

clear all;disp('Clearing');
close all;

disp('Setup');

k=200;

rows=5;
N=rows^2;

grid_rows=5;
grid_cols=5;

N=grid_rows*grid_cols;

B=ones(N,1)*[1,1,1,1,-8,1,1,1,1];
B(1:grid_cols:N,3)=0;%SW
B(grid_cols+1:grid_cols:N,6)=0;%S
B(grid_cols+1:grid_cols:N,9)=0;%SE
B(grid_cols:grid_cols:N,7)=0;%NE
B(grid_cols:grid_cols:N,4)=0;%N
B(grid_cols:grid_cols:N,1)=0;%NW
d=[-grid_cols-1,-grid_cols,-grid_cols+1,-1,0,1,grid_cols-1,grid_cols,grid_cols+1];

A=spdiags(B,d,N,N);

%%
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
fprintf('Iterating:  setup');
for i=1:k
    [xv,WUv]=amg_cycle(A,b,xv,1,7);
    WU=WUv+WU;
    rv=[rv,norm(A*xv-b)];
    fprintf('\b\b\b\b\b\b% 5d\n',i);
end

[xj,rj]=Jacobi(A,b,x,round(WU));

xkn=xv;
rkn=norm(A*xkn-b);
rjn=norm(A*xj-b);

semilogy(linspace(0,length(rj),length(rj)),rj,linspace(0,length(rj),length(rv)),rv);
legend({'Jacobi','AMG cycle'});

