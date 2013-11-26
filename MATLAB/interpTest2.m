function interpTest2
levels=4;
N=127;

Z=peaks(N);
%Z=rand(N);
%Z=Z-mean(mean(Z));
v=reshape(Z',N^2,1); 

interp_grid(levels,v,N);
    function vh = interp_grid(levels,v,nh)
        disp(levels);
        [ Ih_2h , I2h_h ]=get_restriction2D(nh,nh);
        n2h=(nh-1)/2;
        v2h=Ih_2h*v;
        Z2h=reshape(v2h,n2h,n2h)';
        surf(Z2h);
        pause;
        if(levels>0)
            v2h=interp_grid(levels-1,v2h,n2h);
        end
        vh=I2h_h*v2h;
        Zh=reshape(vh,nh,nh)';
        surf(abs(fft2(Zh)));
        pause;
    end
end