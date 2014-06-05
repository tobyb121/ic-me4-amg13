function [ xk , rk] = Jacobi( varargin)
%JACOBI perform iterations of weighted Jacobi
% xk = JACOBI( A, b, x0, k) Perform k iterations of Jacobi on Ax=b
% [ xk , rk] = JACOBI( A, b, x0, k) Perform k iterations of Jacobi on Ax=b
% and return 2-norm of residual with number of iterations (rk)
% JACOBI('w',w) Set weighting
persistent w;

% Check if weight paramter is being set
if(ischar(varargin{1}))
    if strcmp(varargin{1},'w')
        w=varargin{2};
    end
    return;
else
    narginchk(4, 4);
    A=varargin{1};
    b=varargin{2};
    x0=varargin{3};
    k=varargin{4};
end

if isempty(w)
    w=2/3;
end

xk=x0;

% create preconditioner matrix as diagonal part of A divided by w
P=diag(diag(A))/w;

% Invert P (just invert all diagonal elements, but MATLAB does this very
% quickly alreadys so no point looping)
inv_P=P^-1;

% create M=(I-inv_P*A)
M=speye(size(A,1))-inv_P*A;

% if residual is expected in ouput allocate storage
if(nargout==2)
    rk=zeros(k,1);
end

% Do k iterations
for i=1:k
    xk=M*xk+inv_P*b;
    if(nargout==2)
        rk(i)=norm(A*xk-b);
    end
end
end

