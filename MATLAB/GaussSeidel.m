function [ xk , rk] = GaussSeidel( A, b, x0 , k)
xk=x0;
P=diag(diag(A));
inv_P=P^-1;
M=speye(size(A,1))-inv_P*A;

if(nargout==2)
    rk=zeros(k,1);
end

for i=1:k
    for j=1:size(A,1)
        xk(j)=M(j,:)*xk+inv_P(j,:)*b;
    end
    if(nargout==2)
        rk(i)=norm(A*xk-b);
    end
end
end

