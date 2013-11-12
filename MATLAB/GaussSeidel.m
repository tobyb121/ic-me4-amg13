function [ xk ] = GaussSeidel( A, b, x0 , k)
xk=x0;
P=tril(A);
inv_P=P^-1;
M=speye(size(A,1))-inv_P*A;
for i=1:k
    xk=M*xk+inv_P*b;
end
end

