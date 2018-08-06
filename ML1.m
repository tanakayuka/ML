all clear

%最急降下法のステップ幅
lambda =  0.7;
  n = 40;
  omega = randn(1, 1);
  noise = 0.8 * randn(n, 1);
  x = randn(n, 2);
  y = 2 * (omega * x(:, 1) + x(:, 2) + noise > 0) - 1;
  
w=[1,1];
wn=[1,1];
loop = 50;
figure
for i = 1:loop
    psum = zeros(1, 2);
    gradn = zeros(2,1);
    hess = zeros(2,2);
    
       pw = w;
       pwn = wn;
 
    for j = 1:n
        p = 1/(1+exp(-y(j)*(x(j,1)*w(1)+x(j,2)*w(2))));
        psum = psum - (1-p)*y(j)*x(j,:);
        
    end
    
    %最急降下
    grad = psum + lambda*2*w;
    w = w - 0.2*grad;
    
    jw(i) = sum(log(1+exp(-x*pw.'.*y))) + lambda*pw*pw.';
    
    plot(i,jw(i),'bo');
    hold on;
    
end
figure;
for i = 1:loop
    gradn = zeros(2,1);
    hess = zeros(2,2);
       pwn = wn;
 
    for j = 1:n

        pn = 1/(1+exp(-y(j)*(x(j,1)*wn(1)+x(j,2)*wn(2))));    
        gradn = gradn - (1-pn)*y(j)*x(j,:).';
        hess = hess + pn*(1-pn)*x(j,:).'*x(j,:);
    end
    %ニュートン法
    gradn = (gradn) + lambda*2*wn.';
    hess = (hess) + 2*lambda*eye(2);
    d = linsolve(hess,-gradn);
    wn = wn + d.';
    
   
    jwn(i) = sum(log(1+exp(-x*pwn.'.*y))) + lambda*pwn*pwn.';
    
       plot(i,jwn(i),'r*');
       hold on;
       %plot(i,sqrt((w(1)-pw(1))^2+(w(2)-pw(2))^2),'bo');
    
    
end


figure;
for i = 1 :n
    if y(i) == 1
        plot(x(i,1),x(i,2),'ro');
    else
        plot(x(i,1),x(i,2),'bs');
    end
    hold on;
end

z = -3:1/100:3;
plot(z,-wn(1)*z/wn(2),'r')
xlim([-3 3])
ylim([-3 3])

figure;
for i = 1 :n
    if y(i) == 1
        plot(x(i,1),x(i,2),'ro');
    else
        plot(x(i,1),x(i,2),'bs');
    end
    hold on;
end

plot(z,-w(1)*z/w(2),'b')
omega
wn(1)/wn(2)
w(1)/w(2)

%ガウスニュートン法へのリスペクトのやつ