function [ Ih_2h , I2h_h ] = get_restriction2D( grid_rows, grid_cols )

Nh=grid_rows*grid_cols;

rows2h=(grid_rows-1)/2;
cols2h=(grid_cols-1)/2;

N2h=rows2h*cols2h;


n=0;
i=ones(9*N2h,1);
jx=ones(9*N2h,1);
jy=ones(9*N2h,1);
s=zeros(9*N2h,1);

for y=2:2:grid_rows
    for x=2:2:grid_cols
        i(n*9+1:(n+1)*9)=n+1;
        jx(n*9+1:(n+1)*9)=[ x-1, x ,x+1,...
                            x-1, x ,x+1,...
                            x-1, x ,x+1];
        jy(n*9+1:(n+1)*9)=[ y-1,y-1,y-1,...
                             y , y , y ,...
                            y+1,y+1,y+1];
        s(n*9+1:(n+1)*9)=[1,2,1,2,4,2,1,2,1];
        n=n+1;
    end
end
j=(grid_cols*(jy-1))+jx;

Ih_2h=sparse(i,j,s/16,N2h,Nh);
I2h_h=sparse(j,i,s/4,Nh,N2h);
end

