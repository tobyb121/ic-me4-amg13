function [ xh ] = twoGridCorrection( Ah,xh,bh,rows,cols )

v1=3;
v2=10;
v3=3;

%Relax Ax=b v1 times on grid(h)
xh=Jacobi(Ah,bh,xh,2/3,v1);

%Generate restriction operators for grid(h)<=>grid(2h)
[ Ih_2h , I2h_h ]=getRestriction2D(rows,cols);

%Calculate residual
rh=bh-Ah*xh;

%Restrict to grid(2h)
r2h=Ih_2h*rh;
A2h=Ih_2h*Ah*I2h_h;
e2h=zeros(size(r2h));

%Relax Ae=r v2 times on grid(2h)
e2h=Jacobi(A2h,r2h,e2h,2/3,v2);

%Interpolate approximation for e back to grid(h)
eh=I2h_h*e2h;

%Apply error approximation back to x
xh=xh+eh;

%Relax Ax=b v3 times on grid(h)
xh=Jacobi(Ah,bh,xh,2/3,v3);
end

