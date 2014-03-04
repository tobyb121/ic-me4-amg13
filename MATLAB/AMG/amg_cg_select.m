function [C,As]=amg_cg_select(A)
N=size(A,1);

two_pass=1;
theta=0.8;

C=zeros(N,1);
i_c=1;

As=amg_get_strength_matrix(A,theta);

w=sum(As~=0,2)-1;

while(max(w)>0)
    [~,I]=max(w);
    C(i_c)=I;
    i_c=i_c+1;
    new_F=find(As(I,:));
    w(new_F)=0;
    [~,nextC]=find(As(new_F,:));
    for j=nextC'
        if(w(j)>0)
            w(j)=w(j)+1;
        end
    end
end

%% Second Pass
if(two_pass)
    F=1:N;
    F(C(C~=0))=[];
    
    AsFF=As(F,F);
    [rFF,cFF,~]=find(AsFF);
    D=rFF==cFF;
    rFF(D)=[];
    cFF(D)=[];
    
    for i=1:length(rFF)
        s=0;
        for c=C(C~=0)'
            if(As(F(rFF(i)),c)~=0&&As(F(cFF(i)),c)~=0)
                s=1;
                break;
            end
        end
        if(s==0)
            C(i_c)=F(rFF(i));
            i_c=i_c+1;
        end
    end
end

%% De-allocate 0's from C
C(C==0)=[];

end