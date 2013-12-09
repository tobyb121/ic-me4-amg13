function [ xh ] = vcycleDbg(varargin )
%VCYCLE Perform Geometric Multigrid V-Cycle
% xh = VCYCLE(Ah,bh,xh,rows,cols,levels) Solve Ah*xh=bh over 'levels'
% VCYCLE('smoother',func) Set smoother function to use for relaxation steps
% VCYCLE(params...)
% 	'v1',n Set number of pre-restriction smoother steps
% 	'v2',n Set number of smoothing operations on coarsest grid
% 	'v3',n Set number of smoothing steps for post interpolation
persistent smoother v1 v2 v3;

global xk1 xk2 xk3 xk4 rk1 rk2 rk3 rk4;

if(ischar(varargin{1}))
    for i=1:2:nargin
        switch varargin{i}
            case 'v1'
                v1=varargin{i+1};
            case 'v2'
                v2=varargin{i+1};
            case 'v3'
                v3=varargin{i+1};
            case 'smoother'
                smoother=varargin{i+1};
        end
    end
    return;
else
    narginchk(8, 8);
    Ah=varargin{1};
    bh=varargin{2};
    xh=varargin{3};
    rows=varargin{4};
    cols=varargin{5};
    levels=varargin{6};
    xk=varargin{7};
    bk=varargin{8};
end

if(isempty(xk1))
   xk1={};
   xk2={};
   xk3={};
   xk4={};
   rk1={};
   rk2={};
   rk3={};
   rk4={};
end

if isempty(v1)
    v1=10; %pre-restriction
    v2=50; %smooth on coarsest
    v3=10; %post-interpolation
    smoother=@Jacobi;
end

%Relax Ax=b v1 times on grid(h)
xh=smoother(Ah,bh,xh,v1);

xk1{length(xk1)+1}=xk+xh;
rk1{length(rk1)+1}=norm(Ah*xk1{length(xk1)}-bk);
%Generate restriction operators for grid(h)<=>grid(2h)
[ Ih_2h , I2h_h ]=getRestriction2D(rows,cols);
if mod(rows,2) == 0
    rows2h=rows/2;
else
    rows2h=(rows-1)/2;
end 
   
if mod(cols,2) == 0
    cols2h=cols/2;
else
    cols2h=(cols-1)/2;
end
%Calculate residual
rh=bh-Ah*xh;

%Restrict to grid(2h)
r2h=Ih_2h*rh;
A2h=Ih_2h*Ah*I2h_h;
e2h=zeros(size(r2h));
xk2h=Ih_2h*(xk+xh);

xk2{length(xk2)+1}=xk2h;
rk2{length(rk2)+1}=norm(A2h*xk2{length(xk2)}-Ih_2h*bk);
%if not on coarsest grid
if(levels>1)
    %recursive call, to deeper grid
    e2h=vcycleDbg(A2h,r2h,e2h,rows2h,cols2h,levels-1,xk2h,Ih_2h*bk);
else
    %Relax Ae=r v2 times on grid(2h)
    e2h=smoother(A2h,r2h,e2h,v2);
end

xk3{length(xk3)+1}=xk2h+e2h;
rk3{length(rk3)+1}=norm(A2h*xk3{length(xk3)}-Ih_2h*bk);
%Interpolate approximation for e back to grid(h)
eh=I2h_h*e2h;

%Apply error approximation back to x
xh=xh+eh;

xk4{length(xk4)+1}=xk+xh;
rk4{length(rk4)+1}=norm(Ah*xk4{length(xk4)}-bk);

%Relax Ax=b v3 times on grid(h)
xh=smoother(Ah,bh,xh,v3);

end

