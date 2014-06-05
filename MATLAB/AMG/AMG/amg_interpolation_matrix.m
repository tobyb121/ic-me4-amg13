function [Ihc_h]=amg_interpolation_matrix(Ah,C)
%AMG_INTERPOLATION_MATRIX Generates the AMG interpolation matrix
% Ihc_h = amg_interpolation_matrix(Ah, C) Returns the iterpolation matrix
% for A2h=>Ah given C point selection
% amg_interpolation_matrix('reset') Resets the cache

persistent size_cache Ihc_h_cache;

% Check if the cache needs to be reset
if(ischar(Ah))
    if(strcmp(Ah,'reset'))
        Ihc_h_cache={};
        size_cache=[];
    end
    return
end

N=size(Ah,1);

% Check if there is a matrix already in the cache and if so return it
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

% pre-allocate the matrix
Ihc_h=zeros(N,length(C));

% get an array of F points, list all nodes then remove C points
F=1:N;
F(C)=[];

% Start by setting all interpolation weights for C points to 1 (C points
% are interpolated directly)
for c=1:length(C)
    % c is node index on the coarse grid C(c) is node index on fine grid
    Ihc_h(C(c),c)=1;
end

% Calculate interpolation weights for the F points
for f=1:length(F)
    % calculate row sum of coefficients for F-C connections for the F point
    rowsum=sum(abs(Ah(F(f),C)));
    % set the interpolation weights to the absolute value of the
    % coefficient, divided by the rowsum (all F-C connections done at once)
    Ihc_h(F(f),:)=abs(Ah(F(f),C))/rowsum;
end

% convert to sparse storage
Ihc_h=sparse(abs(Ihc_h));

% store in cache
Ihc_h_cache{length(Ihc_h_cache)+1}=Ihc_h;
size_cache(length(size_cache)+1)=size(Ah,1);

end