clear;
N=30;
b=[2,20];
i=linspace(0,1,N);
x=sum(sin(pi*i'*b),2);
xi=sum(sin(pi*linspace(0,1,10*N)'*b),2);

B=[1*ones(N,1),-2*ones(N,1),1*ones(N,1)];
d=[-1,0,1];
A=spdiags(B,d,N,N);

b=A*x;

x0=zeros(N,1);

xn=Jacobi(A,b,x0,2);
plot(x-xn);