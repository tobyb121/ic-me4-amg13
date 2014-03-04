function [r,c]=index2rc(i,N)
   c=ceil(i/N);
   r=mod(i-1,N)+1;
end