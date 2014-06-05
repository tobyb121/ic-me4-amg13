% AMGCYCLEDEPTHTEST AMGCycleDepthTest tests the AMG Cycle against the 80x80
% Poisson problem on cycle depths of 1-7
clear all;
close all;

% Reset caches and initialise V-cycle
amg_cycle('reset');
amg_cycle('v1',3,'v2',10,'v3',3,'smoother',@Jacobi);

% Number of V-cycle iterations
k=500;

% Initialise redisual trace for k iterations on 7 grids
rv=zeros(k,7);

% Load the matrix and RHS, multiply by -1 so that A has positive diagonal,
% negative off diagonal
load('PoissonTestMatrices/poisson_80');
A=-A;
b=-b;

% Initialise initial approximation for x
x=ones(N,1);

for i=1:7
    % Reset the AMG cycle
    amg_cycle('reset');
    
    % Set xv to initial approximation for x
    xv=x;
    
    % Iterate k times
    fprintf('Iterating:  setup');
    for j=1:k
        % Perform AMG V-cycle over i grids
        xv=amg_cycle(A,b,xv,1,i);
        % Calculate the residual and store in the trace
        rv(j,i)=norm(b-A*xv);
        %Print the iteration number
        fprintf('\b\b\b\b\b\b% 5d\n',j);
    end
end

%Plot the results
semilogy(rv);
legend({'1 grid','2 grids','3 grids','4 grids','5 grids','6 grids','7 grids'});