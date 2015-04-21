% 
[m,n]=size(A);
if exist('toler') ~= 1
   toler = 1.e-4;
end;

% set initial variable y
if exist('y0') ~= 1
   y = zeros(m,1);
else
   y=y0;
end;
% set initial slack s
if exist('s0') ~= 1
   s = ones(n,1);
else
   s=s0;
end;
% set initial multiplier x
if exist('x0') ~= 1
   x = ones(n,1);
else
   x=x0;
end;
mu = x'*s/n;
% set ADMM penalty weight beta
if exist('beta') ~= 1
   beta = 1;
end;
% 
[R,p,S] = chol(sparse(A*A' + 10^-7*eye(m)));
assert(p == 0) % is matrix positive definite??? 
%
stuff = 0;
tic
total_it = 0;
progress = [];
gamma = 1 - sqrt(1/n);
while s'*x/sqrt(n) >= toler,
 for k=1:n,
  % Update s
  cc=c-A'*y;
  cc=beta*cc-x;
  s=(sqrt(cc.^2+(4*beta*mu)*ones(n,1))+cc)/(2*beta);
  % Update y
  bb=b+A*(beta*(c-s)-x);
  bb=bb/beta;
  y= S * (R \ (R' \ (S' * bb)));
  % Update multipliers
  x=max(x+beta*(A'*y+s-c),0);
  
  if s'*x/n <= mu*1.5 %(1 + s'*x/(norm(x) + norm(s)))
      disp('break')
      break
  end
  total_it = total_it + 1;
  if total_it > 3*n
      break
  end
  progress(total_it) = s'*x; %/(norm(x)+norm(s));norm(A'*y+s-c)/(1+norm(s))+norm(A*x-b)/(1+norm(x))+
 end;
 if total_it > 3*n
      break
  end
  mu=gamma*mu;
  s'*x/(norm(x) + norm(s))
end;
total_it
toc
norm(A'*y+s-c)/(1+norm(s))+norm(A*x-b)/(1+norm(x))+s'*x


[m,n] = size(A)
% find the actual solution
tic
xTrue = linprog(c,[],[],A,b,zeros(n,1));
toc

tic
xAlmostTrue = linprog(c,[],[],A,b,-ones(n,1));
toc
(c'*xTrue - c'*xAlmostTrue)/(norm(c) + norm(b))
norm(xTrue - xAlmostTrue)/norm(xAlmostTrue)

% primal accuracy
(c'*xTrue - c'*x)/(norm(b) + norm(c))
norm(A*x-b)/(1+norm(x))
norm(min(x,0),2)

% dual accuracy
(b'*y - c'*xTrue)/(norm(b) + norm(y))
norm(A'*y+s-c)/(1+norm(s))

semilogy(progress)
xlabel('number of iterations')
ylabel('Complementary gap =  s^T x')

M = [eye(n) A'; A zeros(m)];
tic
[L,D,p] = ldl(M);
toc

tic
[R,p,S] = chol(A*A');
toc

tic
[R,p,S] = chol(sparse(A*A'+eye(m)));
toc
