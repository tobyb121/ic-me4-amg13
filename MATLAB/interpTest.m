n=31;
Z=peaks(n);
v=reshape(Z',n^2,1);

clf;
subplot(2,3,1);
surf(Z);

[ Ih_2h , I2h_h ]=getRestriction2D(n,n);
[ I2h_4h , I4h_2h ]=getRestriction2D((n-1)/2,(n-1)/2);

v2h=Ih_2h*v;

Z2h=reshape(v2h,(n-1)/2,(n-1)/2)';
subplot(2,3,2);
surf(Z2h);

v4h=I2h_4h*v2h;

Z4h=reshape(v4h,(n+1)/4-1,(n+1)/4-1)';
subplot(2,3,3);
surf(Z4h);

v2h=I4h_2h*v4h;

Z2h=reshape(v2h,(n-1)/2,(n-1)/2)';
subplot(2,3,5);
surf(Z2h);

vh=I2h_h*v2h;
Zh=reshape(vh,n,n)';
subplot(2,3,4);
surf(Zh);