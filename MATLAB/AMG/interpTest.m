clear;

amg_cg_select('reset');
amg_interpolation_matrix('reset');

N=10;

B=0.5*ones(N,1)*[-1,2,-1];
A=spdiags(B,[-1,0,1],N,N);

C=amg_cg_select(A);

Ihc_h=amg_interpolation_matrix(A,C);