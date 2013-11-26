n=127;
Z=peaks(n);
Z=rand(n);
Z=Z-mean(mean(Z));
v=reshape(Z',n^2,1);

h1=figure;
h2=figure;

figure(h1);
subplot(2,3,1);
figure
surf(Z);
set(gca,'visible','off');axis tight;
figure(h2);
subplot(2,3,1);
surf(abs(fft2(Z)));

[ Ih_2h , I2h_h ]=getRestriction2D(n,n);
[ I2h_4h , I4h_2h ]=getRestriction2D((n-1)/2,(n-1)/2);

v2h=Ih_2h*v;

Z2h=reshape(v2h,(n-1)/2,(n-1)/2)';

figure(h1);
subplot(2,3,2);
figure;
surf(Z2h);
set(gca,'visible','off');axis tight;
figure(h2);
subplot(2,3,2);
surf(abs(fft2(Z2h)));

v4h=I2h_4h*v2h;

Z4h=reshape(v4h,(n+1)/4-1,(n+1)/4-1)';

figure(h1);
subplot(2,3,3);
figure;
surf(Z4h);
set(gca,'visible','off');axis tight;
figure(h2);
subplot(2,3,3);
surf(abs(fft2(Z4h)));

v2h=I4h_2h*v4h;

Z2h=reshape(v2h,(n-1)/2,(n-1)/2)';

figure(h1);
subplot(2,3,5);
figure;
surf(Z2h);
set(gca,'visible','off');axis tight;
figure(h2);
subplot(2,3,5);
surf(abs(fft2(Z2h)));

vh=I2h_h*v2h;
Zh=reshape(vh,n,n)';

figure(h1);
subplot(2,3,4);
figure;
surf(Zh);
set(gca,'visible','off');axis tight;
figure(h2);
subplot(2,3,4);
surf(abs(fft2(Zh)));
