clear;
N=8012;
i=1;
B=[-1*ones(N,1),2*ones(N,1),-1*ones(N,1)];
d=[-1,0,1];
A=spdiags(B,d,N,N);
clear('B','d');
x=rand(N,1).*(0.5-0.5*cos((1:N)*pi*2/N))';
b=A*x;
x0=zeros(N,1);
xj=Jacobi(A,b,x0, 1,i);
xjw=Jacobi(A,b,x0, 2/3,i);
xgs=GaussSeidel(A,b,x0,i);
figure;
plot(1:N,x,1:N,xj,1:N,xjw,1:N,xgs);
figure;

e0=x0-x;
ej=xj-x;
ejw=xjw-x;
egs=xgs-x;

s_e0=abs(fft(e0));
s_ej=abs(fft(ej));
s_ejw=abs(fft(ejw));
s_egs=abs(fft(egs));

plot(1:N,s_e0./s_e0,1:N,s_ej./s_e0,1:N,s_ejw./s_e0,1:N,s_egs./s_e0);
legend('e0','Jacobi','Weighted Jacobi','Gauss Seidel');
xlim([3,N/2]);