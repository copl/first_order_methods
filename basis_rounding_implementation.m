[m,n]=size(A);

mubarrierADMM_v2;
%[x,y,s] = ProjectedAGD_V2(c, A, b, n, 0.01 );

xTrue = linprog(c,[],[],A,b,zeros(n,1));

[sorted_x,priority_order] = sort(x,'descend');
priority_order = priority_order(sorted_x > -0.1); %0.1*norm(x)/m);
%priority_order = randperm(n); % random priority order

basis = round_to_basis( A, priority_order );
B = A(:,basis);
%xB = lsqnonneg(B,b);
xB = B \ b;
xExact = zeros(size(A,2),1);
xExact(basis) = xB;



c'*xExact
norm(min(xExact,0))
norm(A*xExact-b)
sum(xExact < 0)

A_prime = A;
A_prime(:,xExact < 0) = -A_prime(:,xExact < 0);
c_prime = c;
c_prime(xExact < 0) = 9999;
c_prime = [c_prime; c(xExact < 0)];
A_prime = [A_prime A(:,xExact < 0)];


bfs = zeros(1,size(A_prime,2));
bfs(basis) = 1;
[xSimplex, iter] = rsm (c_prime', A_prime, b, 0.01, 0.01, 0.01, bfs, 2000);
iter