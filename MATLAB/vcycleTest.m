clear;
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
xj=Jacobi(A,b,x,k);
x2=GaussSeidel(A,b,x,k);
vcycle('v1',10,'v2',50,'v3',10,'smoother',@Jacobi)
xv=vcycle(A,b,x,grid_rows,grid_cols,5);
vcycle('v1',10,'v2',50,'v3',10,'smoother',@GaussSeidel)
xv2=vcycle(A,b,x,grid_rows,grid_cols,5);

Xj=reshape(xj,grid_rows,grid_cols);
X2=reshape(x2,grid_rows,grid_cols);
Xv=reshape(xv,grid_rows,grid_cols);
Xv2=reshape(xv2,grid_rows,grid_cols);

plot_mode=0;
%plot_mode=input('plot mode: ');

if(plot_mode==1)
    ext=input('Extension: *.','s');
    figure('Position',[0,0,900,700]);
    axis tight;caxis([-10,10]);zlim([-10,10]);
    surf_plot(X,'Analytical Solution',ext);
    surf_plot(Xj,'Jacobi (100 iterations)',ext);
    surf_plot(X2,'Two Grid Correction',ext);
    surf_plot(Xv,'Multigrid V-Cycle (5 Grid)',ext);

    surf_plot(Xj-X,['Jacobi (100 iterations) Error (max=',num2str(max(max(abs(Xj-X)))),')'],ext);
    surf_plot(X2-X,['Two Grid Correction Error (max=',num2str(max(max(abs(X2-X)))),')'],ext);
    surf_plot(Xv-X,['Multigrid V-Cycle (5 Grid) Error (max=',num2str(max(max(abs(Xv-X)))),')'],ext);
else
    subplot(2,4,1);surf(X);axis tight
    caxis([-10,10]);zlim([-10,10]);title('Analytical Solution');
    subplot(2,4,2);surf(Xj);axis tight
    caxis([-10,10]);zlim([-10,10]);title('Jacobi (100 iterations)');
    subplot(2,4,3);surf(X2);axis tight
    caxis([-10,10]);zlim([-10,10]);title('Two Grid Correction');
    subplot(2,4,4);surf(Xv);axis tight
    caxis([-10,10]);zlim([-10,10]);title('Multigrid V-Cycle (5 Grid)');

    subplot(2,4,6);surf(Xj-X);axis tight
    caxis([-10,10]);zlim([-10,10]);title(['Jacobi (100 iterations) Error (max=',num2str(max(max(abs(Xj-X)))),')']);
    subplot(2,4,7);surf(X2-X);axis tight
    caxis([-10,10]);zlim([-10,10]);title(['Two Grid Correction Error (max=',num2str(max(max(abs(X2-X)))),')']);
    subplot(2,4,8);surf(Xv-X);axis tight
    caxis([-10,10]);zlim([-10,10]);title(['Multigrid V-Cycle (5 Grid) Error (max=',num2str(max(max(abs(Xv-X)))),')']);
end