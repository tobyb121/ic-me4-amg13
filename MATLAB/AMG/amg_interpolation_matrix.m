function [Ihc_h]=amg_interpolation_matrix(Ah,C)

global Ihc_h_cache; 
persistent size_cache;

if(ischar(Ah))
    if(strcmp(Ah,'reset'))
        Ihc_h_cache={};
        size_cache=[];
    end
    return
end

N=size(Ah,1);

if ~isempty(Ihc_h_cache) 
    for i=1:length(size_cache)
       if(size_cache(i)==N)
           Ihc_h=Ihc_h_cache{i};
          return; 
       end
    end
else
    Ihc_h_cache={};
    size_cache=[];
end

Ihc_h=zeros(N,length(C));

F=1:N;
F(C)=[];

for c=1:length(C)
    Ihc_h(C(c),c)=1;
end

for f=1:length(F)
    rowsum=sum(Ah(F(f),C));
    for c=1:length(C)
        Ihc_h(F(f),c)=Ah(F(f),C(c))/rowsum;
    end
end

Ihc_h=sparse(abs(Ihc_h));

Ihc_h_cache{length(Ihc_h_cache)+1}=Ihc_h;
size_cache(length(size_cache)+1)=size(Ah,1);

end