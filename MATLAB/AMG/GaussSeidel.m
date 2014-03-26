function [ xk , rk] = GaussSeidel( A,b,x0,k)
%% JACOBI perform iterations of weighted Jacobi
% xk = JACOBI( A, b, x0, k) Perform k iterations of Jacobi on Ax=b
% [ xk , rk] = JACOBI( A, b, x0, k) Perform k iterations of Jacobi on Ax=b
% and return 2-norm of residual with number of iterations (rk)
% JACOBI('w',w) Set weighting

xk=x0;
P=diag(diag(A));
inv_P=P^-1;
M=speye(size(A,1))-inv_P*A;

if(nargout==2)
    rk=zeros(k,1);
end

for i=1:k
    for r=1:size(A,1)
        xk(r)=M(r,:)*xk+inv_P(r,:)*b;
    end
    if(nargout==2)
        rk(i)=norm(A*xk-b);
    end
end
end



