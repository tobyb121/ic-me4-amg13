function [ As ] = amg_get_strength_matrix(A,theta)
%AMG_GET_STRENGTH_MATRIX
% As = amg_get_strength_matrix(A,theta) Returns the strength matrix for A
% using a threshold theta: As(i,j)=A(i,j) for strong connections otherwise
% As(i,j)=0. Connection is strong when abs(A(i,j))<theta*max(A:,j)

    % Get list of row, column and values of non-zeros (luckily A is sparse)
    [r,c,v]=find(A);
    
    % Get an array of the maximum coefficient on each row that isn't the
    % diagonal, could potentially just do max of each row, as diagonal
    % should be negative
    aik_max=max(-A+diag(diag(A)));
    
    % Initialise As to A, values can then be removed
    As=A;
    
    % Loop over all non-zeros in the array
    for i=1:length(r)
        if(r(i)~=c(i)) % Ignore diagonal
            
            % If coefficient it less than threshold
            % remove it from As
            if(abs(v(i))<abs(aik_max(r(i))*theta)) 
                As(r(i),c(i))=0;
            end
        end
    end
end

