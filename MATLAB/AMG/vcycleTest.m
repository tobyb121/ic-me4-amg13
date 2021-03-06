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

k=75;

rows=64;
N=rows^2;

grid_rows=rows;
grid_cols=rows;

N=grid_rows*grid_cols;

B=ones(N,1)*[1,1,-4,1,1];
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
xvg=xk0;
rvg=[rk0];
WU=0;
disp('Starting AMG cycle');
fprintf('Iterating:  setup');
for i=1:k
    [xv,WUv]=amg_cycle(A,b,xv,1,7);
    WU=WUv+WU;
    rv=[rv,norm(A*xv-b)];
    xvg=vcycle(A,b,xvg,rows,rows,7);
    rvg=[rvg,norm(A*xvg-b)];
    fprintf('\b\b\b\b\b\b% 5d\n',i);
end

semilogy(1:length(rv),rv,1:length(rvg),rvg);
legend({'AMG Cycle','Geometric V-Cycle'});

