clear all;
clf;

createMGTest;

[C,As]=amg_cg_select(A);

F=setdiff(1:N,C);

[i,j]=index2rc(1:N,rows);

[xf,yf]=index2rc(F,rows);

[x,y]=index2rc(C,rows);


gplot(A,[i',j']);
figure;
hold on;
gplot(As,[i',j']);
plot(x,y,'ro',xf,yf,'k.');
hold off;