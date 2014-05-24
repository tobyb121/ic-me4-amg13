function [x,rcg] = conjgrad(A,b,x,k)
    r=b-A*x;
    p=r;
    rsold=r'*r;
    rcg=zeros(1,k);
    for i=1:k
        Ap=A*p;
        alpha=rsold/(p'*Ap);
        x=x+alpha*p;
        r=r-alpha*Ap;
        rcg(i)=norm(r);
        rsnew=r'*r;
        p=r+rsnew/rsold*p;
        rsold=rsnew;
    end
    rcg(isnan(rcg))=rcg(find(~isnan(rcg),1,'last'));
