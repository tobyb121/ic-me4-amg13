N=8192;
r=rand(N);
r=r-mean(mean(r));
r=r.*((0.5-0.5*cos((1:N)*pi*2/N))'*(0.5-0.5*cos((1:N)*pi*2/N)));
%r2=0.25*(r(1:end-2,1:end-2)+r(1:end-2,3:end)+r(3:end,1:end-2)+r(3:end,3:end));
figure
surf(abs(fft2(r)));
%figure
%surf(abs(fft2(r2)));