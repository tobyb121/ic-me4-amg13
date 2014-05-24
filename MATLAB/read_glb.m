clear;
[f,p]=uigetfile('*.glb');
fID=fopen([p,f]);
l='';

maps=[];
m=1;

while ischar(l)
    l=fgetl(fID);
    if(ischar(l)&&strcmp(strtrim(l),'Mapping array'))
        maps=[maps,m];
    end
    m=m+1;
end
frewind(fID);
disp(['Found: ',num2str(length(maps)),' arrays']);
disp(maps);
start_line=input('Start Line: ');
bandwidth=input('Max Bandwidth: ');
grid_rows=input('Grid Rows: ');
grid_cols=input('Grid Columns: ');
N=grid_rows*grid_cols;

i=zeros(N*bandwidth,1);
j=zeros(N*bandwidth,1);
s=zeros(N*bandwidth,1);
b=zeros(N,1);

for v=1:start_line
    l=fgetl(fID);
end

while ischar(l)&&~strcmp(strtrim(l),'Mapping array')
    l=fgetl(fID);
end

u=1;
for n=1:N
l=fgetl(fID);
v=sscanf(l,'%d',bandwidth+2);
while(length(v)-2<v(2))
    l2=fgetl(fID);
    v=[v;sscanf(l2,'%f',bandwidth+2)];
end
i(u:u+v(2)-1)=v(1);
j(u:u+v(2)-1)=v(3:end);
u=u+v(2);
end

while ischar(l)&&~strcmp(strtrim(l),'Condensed matrix')
    l=fgetl(fID);
end

u=1;
for n=1:N
l=fgetl(fID);
v=sscanf(l,'%f',bandwidth+2);
while(length(v)-2<v(2))
    l2=fgetl(fID);
    v=[v;sscanf(l2,'%f',bandwidth+2)];
end
s(u:u+v(2)-1)=v(3:end);
u=u+v(2);
end

i(u:end)=[];
j(u:end)=[];
s(u:end)=[];

A=sparse(i,j,s,N,N);

while ischar(l)&&~strcmp(strtrim(l),'RHS')
    l=fgetl(fID);
end

for n=1:N
l=fgetl(fID);
v=sscanf(l,'%d %f',2);
b(n)=v(2);
end

fclose(fID);
clear f p fID u v l l2 v n bandwidth start_line i j s m maps;
