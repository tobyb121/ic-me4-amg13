clear;

rows=10;
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