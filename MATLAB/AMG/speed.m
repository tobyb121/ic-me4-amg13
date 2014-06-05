P=diag(diag(A));
for i=1:1000
P=P^-1;
end
for i=1:1000
P=P.^-1;
end