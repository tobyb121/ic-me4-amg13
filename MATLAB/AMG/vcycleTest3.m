clear all;

load('..\\..\\Matrices\\GeometricTest_1_20131203\\matrices.mat');
A=-A4;
b=-b4;

N=length(b);

x=zeros(N,1);

xk0=x;
rk0=norm(A*xk0-b);

k=50;
%%
amg_cycle('v1',2,'v2',10,'v3',2,'smoother',@Jacobi);
xv=xk0;
rv=[rk0];
WU=0;
for i=1:k
    [xv,WUv]=amg_cycle(A,b,xv,1,6);
    WU=WUv+WU;
    rv=[rv,norm(A*xv-b)];
    disp(i);
end

xv2=xk0;
rv2=[rk0];
WU2=0;
amg_cycle('v1',3,'v2',10,'v3',3,'smoother',@GaussSeidel);
while(WU2<WU)
    [xv2,WUv]=amg_cycle(A,b,xv2,1,7);
    WU2=WUv+WU2;
    rv2=[rv2,norm(A*xv2-b)];
end

[xj,rj]=GaussSeidel(A,b,x,round(WU));

[xcg,flag,~,kcg,rcg]=pcg(A,b,1e-14,round(WU));

xkn=xv;



semilogy(1:round(WU),rj,1:round(WU+1),rcg,linspace(0,WU,length(rv)),rv,linspace(0,WU,length(rv2)),rv2);
legend({'Jacobi','CG','AMG1','AMG2'});
