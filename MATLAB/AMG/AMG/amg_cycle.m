function [ xh, WU ] = amg_cycle(varargin)
%AMG_CYCLE Perform Algebraic Multigrid Cycle
% xh = amg_cycle(Ah,bh,xh,mu,maxlevels) Solve Ah*xh=bh over 'maxlevels'
% amg_cycle('smoother',func) Set smoother function to use for relaxation steps
% amg_cycle(params...)
% 	'v1',n Set number of pre-restriction smoother steps
% 	'v2',n Set number of smoothing operations on coarsest grid
% 	'v3',n Set number of smoothing steps for post interpolation
% amg_cycle('reset') Reset cache
persistent smoother v1 v2 v3;
persistent size_cache A2h_cache;

% Process Input arguments
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
            case 'reset' % Reset caches
                amg_cg_select('reset');
                amg_interpolation_matrix('reset');
                A2h_cache={};
                size_cache=[];      
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

WU=0;

%Relax Ax=b v1 times on grid(h)
xh=smoother(Ah,bh,xh,v1);

WU=WU+v1;


C=amg_cg_select(Ah);

lenC=length(C);

%Generate interpolation operator for grid(hc)=>grid(h)
[ Ihc_h ]=amg_interpolation_matrix(Ah,C);

% if Ihc_c is empty cannot coarsen further
if isempty(Ihc_h)
    %Relax Ax=b v2 times and return
    xh=smoother(Ah,bh,xh,v2);
    WU=WU+v2;
    return;
end

max_row_Ih_hc=max(sum(Ihc_h,1)); % calculate max row sum

% create Restriction operator for grid(h)=>grid(hc) using variational
% property, scale by 1/max row sum to maintain magnitude between grids
Ih_hc=1/max_row_Ih_hc*Ihc_h';


A2h=[];
% try an load A2h from cache
if ~isempty(A2h_cache)
    for i=1:length(size_cache)
       if(size_cache(i)==lenC)
           A2h=A2h_cache{i};
       end
    end
else
    A2h_cache={};
    size_cache=[];
end
% if it wasn't loaded not in cache so must be generated
if (isempty(A2h))
    A2h=Ih_hc*Ah*Ihc_h;
    A2h_cache{length(A2h_cache)+1}=A2h; % added to cache
    size_cache(length(size_cache)+1)=lenC;
end

for m=1:mu
    %Calculate residual
    rh=bh-Ah*xh;
    
    %Restrict to grid(2h)
    r2h=Ih_hc*rh;
    e2h=zeros(size(r2h));

    %if not on coarsest grid
    if(maxlevels>1)
        %recursive call, to deeper grid
        [e2h,WUc]=amg_cycle(A2h,r2h,e2h,mu,maxlevels-1);
    else
        %Relax Ae=r v2 times on grid(2h)
        e2h=smoother(A2h,r2h,e2h,v2);
        WUc=v2;
    end
    
    %scale work units from the coarse grid to the fine grid
    WU=WU+WUc*length(e2h)/length(xh);
    
    %Interpolate approximation for e back to grid(h)
    eh=Ihc_h*e2h;

    %Apply error approximation back to x
    xh=xh+eh;

    %Relax Ax=b v3 times on grid(h)
    xh=smoother(Ah,bh,xh,v3);
    WU=WU+v3;
end

end

