% 
[m,n]=size(A);
if exist('toler') ~= 1
   toler = 1.e-5;
end;
if exist('gamma') ~= 1
   gamma = 0.9;
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
mu =x'*s/n;
% set ADMM penalty weight beta
if exist('beta') ~= 1
   beta = 1;
end;
% 
[R,p,S] = chol(sparse(A*A' + 10^-5*eye(m)));
assert(p == 0) % is matrix positive definite??? 
%
total_it = 0;
while mu >= toler,
 for k=1:sqrt(1/mu),
  % Update y
  bb=b+A*(beta*(c-s)-x);
  bb=bb/beta;
  y= S * (R \ (R' \ (S' * bb)));
  % Update s
  cc=c-A'*y;
  cc=beta*cc-x;
  s=(sqrt(cc.^2+(4*beta*mu)*ones(n,1))+cc)/(2*beta);
  % Update multipliers
  x=x+beta*(A'*y+s-c);
  total_it = total_it + 1;
 end;
  mu=gamma*mu

end;
total_it
norm(A'*y+s-c)/(1+norm(s))+norm(A*x-b)/(1+norm(x))+s'*x/(norm(x)+norm(s))

% find the actual solution
tic
xTrue = linprog(c,[],[],A,b,zeros(n,1));
toc

% primc'*xTrue - c'*xal accuracy
norm(A*x-b)/(1+norm(x))

% dual accuracu
(b'*y - c'*xTrue)/norm(x)
norm(A'*y+s-c)/(1+norm(s))


% ADMM for solving the barrier problem
%
%      minimize      -b'*y-mu B(s)
%      subject to     A'*y + s = c
%
%      Input: A, b, c, mu
%             and any initial x0>0 and s0 > 0, and multiplier y0.
%      Output: (y,s) of the dual, and x for the primal 
%          
%