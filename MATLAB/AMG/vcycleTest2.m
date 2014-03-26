%% vcycleTest2
% Performs iterations of amg v-cycle on a 64x64 element grid on discretised
% equation:
% 
% a*d2f/dx^2+b*d2f/dy^2+c*d2f/dxdy=0
%        
%           |
%   a=1     |  a=1
%   b=1000  |  b=1
%   c=0     |  c=2
% __________|_________
%           |
%   a=1     |  a=1000
%   b=1     |  b=1
%   c=0     |  c=0
%           |
% Compares result to iteration of Jacobi (for same number of WU) and plots
% residuals for both cases

clear all;
close all;

k=50;

rows=64;
N=rows^2;

n=@(i,j)i+(j-1)*rows;

Bxx=ones(N,1)*[1,-2,1];
Byy=ones(N,1)*[1,-2,1];
Bxy=ones(N,1)*[1,1,-4,1,1];

Bxx(rows:rows:N,1)=0;
Bxx(rows+1:rows:N,3)=0;

Bxy(1:rows:N,2)=0;
Bxy(1:rows:N,5)=0;

dxx=[-1,0,1];
dyy=[-rows,0,rows];
dxy=[-rows-1,-rows+1,0,rows-1,rows+1];

a=zeros(N,1);
b=zeros(N,1);
c=zeros(N,1);
for i=ceil(1:rows/2)
 for j=ceil(1:rows/2)
     a(n(i,j))=-1;
     b(n(i,j))=-1;
     c(n(i,j))=0;
 end
end
for i=ceil(rows/2+1:rows)
 for j=ceil(1:rows/2)
     a(n(i,j))=-1000;
     b(n(i,j))=-1;
     c(n(i,j))=0;
 end
end
for i=ceil(1:rows/2)
 for j=ceil(rows/2+1:rows)
     a(n(i,j))=-1;
     b(n(i,j))=-1000;
     c(n(i,j))=0;
 end
end
for i=ceil(rows/2+1:rows)
 for j=ceil(rows/2+1:rows)
     a(n(i,j))=-1;
     b(n(i,j))=-1;
     c(n(i,j))=-2;
 end
end

Axx=spdiags(a,0,N,N)*spdiags(Bxx,dxx,N,N);
Ayy=spdiags(b,0,N,N)*spdiags(Byy,dyy,N,N);
Axy=spdiags(c,0,N,N)*spdiags(Bxy,dxy,N,N);

A=Axx+Ayy+Axy;

clear(  'Bxx','Byy','Bxy',...
        'dxx','dyy','dxy',...
        'Axx','Ayy','Axy',...
        'a','b','c');

u=zeros(N,1);

b=A*u;
x=rand(N,1);

xk0=x;
rk0=norm(A*xk0-b);
%%
amg_cycle('v1',3,'v2',10,'v3',3,'smoother',@Jacobi);
xv=xk0;
rv=[rk0];
WU=0;
fprintf('Iterating:  setup');
for i=1:k
    [xv,WUv]=amg_cycle(A,b,xv,1,6);
    WU=WUv+WU;
    rv=[rv,norm(A*xv-b)];
    fprintf('\b\b\b\b\b\b% 5d\n',i);
end

[xj,rj]=Jacobi(A,b,x,round(WU));

xkn=xv;
rkn=norm(A*xkn-b);
rjn=norm(A*xj-b);

semilogy(linspace(0,length(rj),length(rj)),rj,linspace(0,length(rj),length(rv)),rv);

