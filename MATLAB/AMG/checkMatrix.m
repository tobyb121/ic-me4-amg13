createMGTest;

A=full(A);

X=reshape(1:N,sqrt(N),sqrt(N));

for x=1:N
   [r,c]=find(X==x);
   [r2,c2]=index2rc(x,sqrt(N));
   disp([r,c,r2,c2]);
end

pause

for i=1:N
   disp(reshape(A(i,:),sqrt(N),sqrt(N)));
   pause;
   
end