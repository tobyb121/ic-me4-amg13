function [ Ih_2h ] = get_restriction2D( grid_rows, grid_cols )

Nh=grid_rows*grid_cols;

rows2h=(grid_rows-1)/2+1;
cols2h=(grid_cols-1)/2+1;

N2h=rows2h*cols2h;


k=1;
bi=ones(N2h*9,1);
bj=ones(N2h*9,1);
bs=zeros(N2h*9,1);

% for i=1:2:grid_rows
%     for j=1:2:grid_cols
%         i2h=(i-1)/2+1;
%         j2h=(j-1)/2+1;
%         v=(i-1)*grid_cols+j;
%         v2h=(i2h-1)*cols2h+j2h;
%         bi(k:k+8)=v2h;
%         bj(k:k+8)=[[v-1, v, v+1]-grid_cols,...
%                     v-1, v, v+1,...
%                    [v-1, v, v+1]+grid_cols];
%         bs(k:k+8)=[ 1,2,1,...
%                     2,4,2,...
%                     1,2,1];
%         k=k+9;
%     end
% end

for i=2:2:grid_rows+1
    for j=2:2:grid_cols+1
        i2h=i/2;
        j2h=j/2;
        v=(i-1)*(grid_cols+2)+j;
        v2h=(i2h-1)*cols2h+j2h;
        bi(k:k+8)=v2h;
        bj(k:k+8)=[[v-1, v, v+1]-grid_cols-2,...
                    v-1, v, v+1,...
                   [v-1, v, v+1]+grid_cols+2];
        bs(k:k+8)=[ 1,2,1,...
                    2,4,2,...
                    1,2,1];
        k=k+9;
    end
end

%I=(bi<=0)|(bi>N2h)|(bj<=0)|(bj>Nh);
%bi(I)=[];
%bj(I)=[];
%bs(I)=[];

Ih_2h=sparse(bi,bj,bs);



end

