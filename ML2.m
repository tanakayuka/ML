clear all;

A = [3 0.5;0.5 1];
eta = 1/(2*max(eig(A)));
w_hat = [0.82;1.09]; %lambda = 2
%w_hat = [0.64;0.18]; %lambda = 4
%w_hat = [0.33;0.]; %lambda = 6
lambda = 2;
wpre = [3;-1];
wnex = [0;0];
mu = [1;2];
par_mu = [0;0];
q = eta*lambda;
ans = zeros(100,2);
for i=1: 100
    semilogy(i,norm(wpre - w_hat),'bo');
    %plot(wpre(1),wpre(2),'bo');
    hold on
    grad = 2*A*(wpre-mu)
    par_mu = wpre - eta*grad;
    for j=1:2
        if par_mu(j) > q
            wnex(j) = par_mu(j) - q;
        elseif par_mu(j) < -q
            wnex(j) = par_mu(j) + q;
        else
            wnex(j) = 0;
        end
    end
    ans(i,:) = wpre;
    wpre = wnex;
end
figure;
plot(ans(:,1),ans(:,2),'bo-')
xlim([-1.5 3])
ylim([-1.5 3])