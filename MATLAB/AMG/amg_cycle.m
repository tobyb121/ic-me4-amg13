function [ xh ] = amg_cycle(varargin)
%AMG_AMG_CYCLE Perform Algebraic Multigrid Cycle
% xh = AMG_CYCLE(Ah,bh,xh,mu,maxlevels) Solve Ah*xh=bh over 'maxlevels'
% AMG_CYCLE('smoother',func) Set smoother function to use for relaxation steps
% AMG_CYCLE(params...)
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
    narginchk(5, 5);
    Ah=varargin{1};
    bh=varargin{2};
    xh=varargin{3};
    mu=varargin{4};
    maxlevels=varargin{5};
end

if isempty(v1)
    v1=3; %pre-restriction
    v2=10; %smooth on coarsest
    v3=3; %post-interpolation
    smoother=@Jacobi;
end

%Relax Ax=b v1 times on grid(h)
xh=smoother(Ah,bh,xh,v1);

C=amg_cg_select(Ah);

%Generate restriction operators for grid(h)<=>grid(hc)
[ Ihc_h ]=amg_interpolation_matrix(Ah,C);

Ih_hc=Ihc_h';

Ih_hcs=sum(Ih_hc,2);

for n=1:length(Ih_hcs)
   Ih_hc(n,:)=Ih_hc(n,:)*1/Ih_hcs(n);
end

%Calculate residual
rh=bh-Ah*xh;

for m=1:mu
    %Restrict to grid(2h)
    r2h=Ih_hc*rh;
    A2h=Ih_hc*Ah*Ihc_h;
    e2h=zeros(size(r2h));

    %if not on coarsest grid
    if(maxlevels>1)
        %recursive call, to deeper grid
        e2h=amg_cycle(A2h,r2h,e2h,mu,maxlevels-1);
    else
        %Relax Ae=r v2 times on grid(2h)
        e2h=smoother(A2h,r2h,e2h,v2);
    end

    %Interpolate approximation for e back to grid(h)
    eh=Ihc_h*e2h;

    %Apply error approximation back to x
    xh=xh+eh;

    %Relax Ax=b v3 times on grid(h)
    xh=smoother(Ah,bh,xh,v3);
end

end

