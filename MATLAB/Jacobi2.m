function [ xk , rk] = Jacobi2( varargin)
%JACOBI perform iterations of weighted Jacobi
% xk = JACOBI( A, b, x0, k) Perform k iterations of Jacobi on Ax=b
% [ xk , rk] = JACOBI( A, b, x0, k) Perform k iterations of Jacobi on Ax=b
% and return 2-norm of residual with number of iterations (rk)
% JACOBI('w',w) Set weighting
persistent w;

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
P=diag(diag(A))/w;
inv_P=P^-1;
M=speye(size(A,1))-inv_P*A;

rk=zeros(k,1);

for i=1:k
    xk=M*xk+inv_P*b;
    rk(i)=norm(A*xk-b);
end
plot(rk);
end

