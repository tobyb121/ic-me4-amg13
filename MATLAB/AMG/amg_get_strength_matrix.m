function [ As ] = amg_get_strength_matrix(A,theta)
    [r,c,v]=find(A);
    
    aik_max=max(-A+diag(diag(A)));
    
    As=A;
    
    for i=1:length(r)
        if(r(i)~=c(i))
            if(-v(i)<aik_max(r(i))*theta)
                As(r(i),c(i))=0;
            end
        end
    end

%     for i=1:size(A,1)
%         for k=(c(r==i&c~=i))';
%             if -A(i,k)<theta*aik_max(i)
%                 As(i,k)=0;
%             end
%         end
%     end
end

