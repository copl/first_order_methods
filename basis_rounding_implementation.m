[m,n]=size(A);
xTrue = linprog(c,[],[],A,b,zeros(n,1));

[sorted_x,priority_order] = sort(x,'descend' );
%priority_order = priority_order(sorted_x > -0.1); %0.1*norm(x)/m);

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



rsm (c, A, b, 0.01, 0.01, 0.01, basis)