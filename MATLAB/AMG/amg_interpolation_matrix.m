function [Ih_hc]=amg_interpolation_matrix(Ah,C)

persistent Ih_hc_cache size_cache;

if ~isempty(Ih_hc_cache) 
    for i=length(size_cache)
       if(s(i)==size(Ah,1))
          Ih_hc=Ih_hc_cache(i);
          return; 
       end
    end
else
    Ih_hc={};
end



Ih_hc_cache{length(Ih_hc_cache)+1}=Ih_hc;
size_cache(length(size_cache)+1)=size(Ah,1);

end