n = 40;
omega = randn(1, 1);
noise = 0.8 * randn(n, 1);
x = randn(n, 2);
y = 2 * (omega * x(:, 1) + x(:, 2) + noise > 0) - 1;
eta = 0.01;
K = zeros(40,40);
for i = 1 : 40
for j = 1:40
K(i,j) = y(i)*y(j)*x(i,:)*x(j,:).';
end
end
lambda = 1;
alpha = ones(40,1);
new_alpha = zeros(40,1);
figure;
for i=1:100
new_alpha = alpha - eta*(((1/(2*lambda))*K*alpha)-ones(40,1));
for j=1:40
if new_alpha(j) < 0
new_alpha(j) = 0;
elseif new_alpha(j) > 1
new_alpha(j) = 1;
end
end
alpha = new_alpha;
w = zeros(2,1);
for j=1 : 40
w = w + alpha(j)*y(j)*x(j,:).';
end
sum = 0;
for j=1:40
if 1 - y(j)*(x(j,:)*w) > 0
sum = sum + 1 - y(j)*(x(j,:)*w);
end
end
sum = sum + lambda* (w.'*w);
plot(i,sum,'bo')
hold on;
plot(i,(-1/(4*lambda))*alpha.'*K*alpha + ones(1,40)*alpha ,'r*');
end
legend('sum of hinge loss function','score of the dual Lagrange function')
figure;
for i=1:n
if y(i) == 1
plot(x(i,1),x(i,2),'bo')
else
plot(x(i,1),x(i,2),'r*')
end
hold on
end
z = -3:1/10:3;
plot(z,-(w(1)/w(2))*z)
