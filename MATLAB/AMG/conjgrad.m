function [xk,rcg] = conjgrad(A,b,x,k)
%CONJGRAD Iterative Conjugate Gradient Solver
% [xk,rcg] = conjgrad(A,b,x,k) Performs k iterations of the conjugate
% gradient method on Ax=b, returning new  approximation for x and the norm
% of the residual at each iteration
% For details of algorithm see:
% http://en.wikipedia.org/wiki/Conjugate_gradient_method#The_resulting_algorithms
    xk=x;
    % calculate residual
    r=b-A*xk;
    p=r;
    rcg=zeros(1,k);
    for i=1:k
        Ap=A*p;
        alpha=(r'*r)/(p'*Ap);
        xk=xk+alpha*p;
        r_new=r-alpha*Ap;
        rcg(i)=norm(r_new);
        beta=(r_new'*r_new)/(r_old'*r_old);
        p=r_new+beta*p;
        r=r_new;
    end
    rcg(isnan(rcg))=rcg(find(~isnan(rcg),1,'last'));
