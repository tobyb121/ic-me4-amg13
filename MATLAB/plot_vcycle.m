global xk1 xk2 xk3 xk4 rk1 rk2 rk3 rk4;

levels=3;

for i=1:levels
    subplot(length(xk1),levels,levels*(i-1)+1);
    XK=reshape(xk1{i},sqrt(size(xk1{i},1)),sqrt(size(xk1{i},1)));
    surf(XK);   
    axis tight;
    caxis([-10,10]);
    %zlim([-10,10]);
    %set(gca,'visible','off');
    title(num2str(rk1{i}));
end

for i=1:levels
    subplot(length(xk2),4,levels*(i-1)+2);
    XK=reshape(xk2{i},sqrt(size(xk2{i},1)),sqrt(size(xk2{i},1)));
    surf(XK);   
    axis tight;
    caxis([-10,10]);
    %zlim([-10,10]);
    %set(gca,'visible','off');
    title(num2str(rk2{i}));
end

for i=1:levels
    subplot(length(xk3),4,levels*(i-1)+3);
    XK=reshape(xk3{i},sqrt(size(xk3{i},1)),sqrt(size(xk3{i},1)));
    surf(XK);   
    axis tight;
    caxis([-10,10]);
    %zlim([-10,10]);
    %set(gca,'visible','off');
    title(num2str(rk3{i}));
end

for i=1:levels
    subplot(length(xk4),4,levels*(i-1)+4);
    XK=reshape(xk4{i},sqrt(size(xk4{i},1)),sqrt(size(xk4{i},1)));
    surf(XK);   
    axis tight;
    caxis([-10,10]);
    %zlim([-10,10]);
    %set(gca,'visible','off');
    title(num2str(rk4{i}));
end
