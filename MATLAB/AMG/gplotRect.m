function gplotRect( A, rows )
[i,j]=index2rc(1:rows^2,rows);
gplot(A,[i',j']);
end

