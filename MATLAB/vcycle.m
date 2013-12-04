function [ xh ] = vcycle(varargin )
%VCYCLE Perform Geometric Multigrid V-Cycle
% xh = VCYCLE(Ah,bh,xh,rows,cols,levels) Solve Ah*xh=bh over 'levels'
% VCYCLE('smoother',func) Set smoother function to use for relaxation steps
% VCYCLE(params...)
% 	'v1',n Set number of pre-restriction smoother steps
% 	'v2',n Set number of smoothing operations on coarsest grid
% 	'v3',n Set number of smoothing steps for post interpolation
persistent smoother v1 v2 v3;

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
    narginchk(6, 6);
    Ah=varargin{1};
    bh=varargin{2};
    xh=varargin{3};
    rows=varargin{4};
    cols=varargin{5};
    levels=varargin{6};
end

disp(rows);

if isempty(v1)
    v1=10; %pre-restriction
    v2=50; %smooth on coarsest
    v3=10; %post-interpolation
    smoother=@Jacobi;
end

%Relax Ax=b v1 times on grid(h)
xh=smoother(Ah,bh,xh,v1);

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

%if not on coarsest grid
if(levels>1)
    %recursive call, to deeper grid
    e2h=vcycle(A2h,e2h,r2h,rows2h,cols2h,levels-1);
else
    %Relax Ae=r v2 times on grid(2h)
    e2h=smoother(A2h,r2h,e2h,v2);
end

%Interpolate approximation for e back to grid(h)
eh=I2h_h*e2h;

%Apply error approximation back to x
xh=xh+eh;

%Relax Ax=b v3 times on grid(h)
xh=smoother(Ah,bh,xh,v3);

end

