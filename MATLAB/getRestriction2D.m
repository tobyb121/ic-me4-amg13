function [ Ih_2h , I2h_h ] = getRestriction2D( grid_rows, grid_cols )

Nh=grid_rows*grid_cols;

if mod(grid_rows,2) == 0
    rows2h=grid_rows/2;
else
    rows2h=(grid_rows-1)/2;
end 
   
if mod(grid_cols,2) == 0
    cols2h=grid_cols/2;
else
    cols2h=(grid_cols-1)/2;
end

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

outside=jx>grid_cols|jy>grid_rows;
i(outside)=[];
jx(outside)=[];
jy(outside)=[];
s(outside)=[];

j=(grid_cols*(jy-1))+jx;

Ih_2h=sparse(i,j,s/16,N2h,Nh);
I2h_h=sparse(j,i,s/4,Nh,N2h);
end

