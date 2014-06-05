% Compares the AMG solver against a Jacobi and Conjugate Gradient solver
% for the Poisson test problem for square grids with 20,40,50,60,70&80
% nodes in each direction. 75 AMG cycles are performed and the equivalent
% number of Jacobi/CG iterations is performed (based on Work Units)

clear all;
close all;

% Grid sizes to load
G=[20,40,50,60,70,80];

% V-cycle iterations
k=75;

% Pre-allocate residual traces for k iterations on each of the G grids
rv=zeros(k,length(G));
rj=zeros(k,length(G));
rcg=zeros(k,length(G));

for i=1:length(G)
    % Load the matrix from the file
    load(['PoissonTestMatrices/poisson_',num2str(G(i))]);
    % A should have negative diagonal and positive off-diagonal elements so
    % needs to be multiplied by -1, same for b so that x has correct sign
    A=-A;
    b=-b;
    % Initial approximation for x
    x=ones(N,1);
       
    % Reset caches and initialise V-cycle
    amg_cycle('reset');
    amg_cycle('v1',3,'v2',10,'v3',3,'smoother',@Jacobi);
    
    % AMG
    %Initialise AMG to initial approximation
    xv=x;
    WU=0;
    fprintf('Iterating:  setup');
    for j=1:k
        [xv,WUv]=amg_cycle(A,b,xv,1,5); % Perform AMG V-cycle iteration over 5 levels
        
        WU=WUv+WU; % Increment total work done by V-cycles
        rv(j,i)=norm(b-A*xv); % calculate residual for this iteration
        fprintf('\b\b\b\b\b\b% 5d\n',j); % Print iteration number
    end
    
    %Jacobi
    fprintf('Iterating: Jacobi %d\n',WU);
    [xj,rjn]=Jacobi(A,b,x,k*ceil(WU/k)); % Iterate Jacobi for same number of work units
    rjn=rjn(1:length(rjn)/k:length(rjn));% Rescale iterations to equivalent V-cycle iterations
    rj(:,i)=rjn;% store residuals for this grid
    
    %Conjugate Gradient
    fprintf('Iterating: CG %d\n',WU);
    [xcg,rcgn]=conjgrad(A,b,x,k*ceil(WU/k)); % Iterate Conjugate Gradient for same number of work units
    rcgn=rcgn(1:length(rcgn)/k:length(rcgn));% Rescale iterations to equivalent V-cycle iterations
    rcg(:,i)=rcgn;% store residuals for this grid
    
end

% Plot results
for g=1:length(G)
    subplot(2,3,g);
    semilogy([rv(:,g),rj(:,g),rcg(:,g)]);
    title(['Grid: ',num2str(G(g))]);
    legend({'AMG','Jacobi','Conjugate Gradient'});
end