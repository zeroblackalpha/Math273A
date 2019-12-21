rng(0);
n = 5;
epsilon = 1e-6;
% A = eye(n);
A = rand(n);
b = 10*rand(n,1);
c = -10*rand(n,1);

tic
[x_simplex,fval,exitflag,output]=linprog(c,A,b,[],[],zeros(n,1),b);
output
toc

% tic
% L = 100;
% x0 = rand(n,1);
% x = zeros(n,1);
% for i=1:n
%     if x0(i) > b(i)
%         x(i) = b(i);
%     else
%         x(i) = x0(i);
%     end
% end
% iter = 0;
% while norm(x-b) > epsilon
%     iter = iter + 1;
%     x_hat = x - L*c;
%     for i=1:n
%         if x_hat(i) > b(i)
%             x(i) = b(i);
%         elseif x_hat(i) < 0
%             x(i) = 0;
%         else
%             x(i) = x_hat(i);
%         end
%     end
% end
% iter
% toc

tic
sigma = sqrt(1/(norm(A)*norm(A)));
tau = sigma;
theta = 1;
x_new = rand(n,1);
y_new = rand(n,1);
iter = 0;
% while norm(x_new-b) > epsilon
x_old = x_new;
y_old = y_new;
iter = iter + 1;
x_hat = x_old- tau*transpose(A)*y_old;
x_new = max(0, x_hat-tau*c);
x_bar = x_new + theta*(x_new-x_old);
y_hat = y_old + sigma*A*x_bar;
y_new = max(0, y_hat-sigma*b);
while norm(x_new-x_simplex) > epsilon
    x_old = x_new;
    y_old = y_new;
    iter = iter + 1;
    x_hat = x_old- tau*transpose(A)*y_old;
    x_new = max(0, x_hat-tau*c);
    x_bar = x_new + theta*(x_new-x_old);
    y_hat = y_old + sigma*A*x_bar;
    y_new = max(0, y_hat-sigma*b);
end
iter
toc

options = optimoptions('linprog','Algorithm','interior-point');
tic
[x,fval,exitflag,output]=linprog(c,A,b,[],[],zeros(n,1),b,options);
output
toc