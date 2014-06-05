function [C]=amg_cg_select(A)
%AMG_CG_SELECT Perform algebric multigrid coarse grid selection
% [C] = amg_cg_select(A) Generate coarse grid selection for matrix A,
% returns list of node indices
% amg_cg_select('reset') Reset cache
persistent size_cache C_cache;

% Check if cache should be cleared
if(ischar(A))
    if(strcmp(A,'reset'))
        C_cache={};
        size_cache=[];
    end
    return
end

% check if there is a selection in the cache already
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

% Use two pass algorithm? (1/0)
two_pass=1;

theta=0.4;

% Pre-allocate space for C-point selection
C=zeros(N,1);
i_c=1; % indexer used for position in C array
 
% Get strength matrix
As=amg_get_strength_matrix(A,theta);

% Get an array of the number the of non-zeros on each column (influence
% value), the minus one is to remove the influence of the node on itself
% (caused by the diagonal)
w=sum(As~=0,2)-1;

%while there are still nodes to select
while(max(w)>0)
    [~,I]=max(w); % get the index of the node with max influence
    C(i_c)=I; % add this node to the C points
    i_c=i_c+1;
    new_F=find(As(I,:)); % Add all connected nodes to F
    w(new_F)=0; % remove the new F nodes from the search list
    
    % Get all nodes connected to newly selected F points as candidates for
    % next C point selection, next C is an array of the column indexes for
    % all non-zeros in As the rows of new_F, ie list of all nodes strongly
    % connected to newly added F points
    [~,nextC]=find(As(new_F,:));
    
    for j=nextC'
        if(w(j)>0) % If node hasn't already been added to C or F
            % Increment influence so it is more likely to be selected on
            % next iteration
            w(j)=w(j)+1;
        end
    end
end

%% Second Pass
if(two_pass)
    % Get list of F points
    F=1:N; % Assume all points are F points
    F(C(C~=0))=[]; % Remove any points that are C points
    
    % Get sub-matrix of As, for F points only (only select rows and columns
    % of the F points)
    AsFF=As(F,F);
    
    % Get the row and column indices of all F-F connections (non-zeros
    % elements of AsFF)
    [rFF,cFF,~]=find(AsFF);
    % Remove diagonal elements (F points connected to themselves)
    D=rFF==cFF;
    rFF(D)=[];
    cFF(D)=[];
    
    % Loop over all F-F connections
    for i=1:length(rFF)
        if(rFF(i)~=0) % If this point was already made into a C point
            
            c=nonzeros(C); % get non-zero elements of C
            
            % If neither of the F points share a connection (gets a list of
            % coefficients  in columns of current C points for the rows two
            % F points and checks if any non-zeros in both rows)
            if(~any(As(F(rFF(i)),c)&As(F(cFF(i)),c)))
                % If there is no shared F-F connections add this point to C
                C(i_c)=F(rFF(i));
                % no need to check other F-F connections now that it's a C
                % point
                rFF(rFF==rFF(i))=0;
                i_c=i_c+1;
            end
        end
    end
end

% De-allocate 0's from C
C(C==0)=[];
C=sort(C);

% Add C point selection to cache
C_cache{length(C_cache)+1}=C;
size_cache(length(size_cache)+1)=size(A,1);
end