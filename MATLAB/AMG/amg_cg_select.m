function [C,As]=amg_cg_select(A)

global C_cache;
persistent size_cache;

if(ischar(A))
    if(strcmp(A,'reset'))
        C_cache={};
        size_cache=[];
    end
    return
end

if ~isempty(C_cache)
    for i=1:length(size_cache)
        if(size_cache(i)==size(A,1))
            C=C_cache{i};
            return;
        end
    end
else
    C_cache={};
    size_cache=[];
end

N=size(A,1);

two_pass=1;
theta=0.4;

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
        if(rFF(i)~=0)
            c=nonzeros(C);
            if(~any(As(F(rFF(i)),c)&As(F(cFF(i)),c)))
                C(i_c)=F(rFF(i));
                rFF(rFF==rFF(i))=0;
                i_c=i_c+1;
            end
        end
    end
end

%% De-allocate 0's from C
C(C==0)=[];
C=sort(C);

C_cache{length(C_cache)+1}=C;
size_cache(length(size_cache)+1)=size(A,1);
end