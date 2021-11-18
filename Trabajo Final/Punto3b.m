c = 20;
m = 20;
f = @(x)exp(-x);
g = @(x)-2*exp(-x);
l = @(t)exp(-2*t);
r = @(t) exp(-1-2*t);
h = 1/m;
k = h/(c*2); 
sigma = c*k/h;
n = m;
a = 0;
b =1;
at = 0;
bt = 1;
X = a+(0:m)*h; %Vector de xs
T = b+(0:n)*k; %Vector de ts
N =m+1;
M = n+1;
A = diag((2-2*(sigma^2))*ones(1,M)) + diag((sigma^2)*ones(1,M-1),-1)+ diag((sigma^2)*ones(1,M-1),1); %Matriz tridiagonal
Fx = ones(M,1);
Gx = ones(M,1);
T0 = zeros(M,1);
T0(1) = l(0);
T0(M) = r(0);

for i=1:M
   Fx(i,1) = f(X(i)); 
   Gx(i,1) = g(X(i));
end

W0 = Fx;
W1 = (1/2)*A*Fx +k*Gx+(sigma^2)/2*T0;
Wij = ones(M,N);
Wij(:,1) = W0;
Wij(:,2) = W1;
for i=3:N
    for j=2:M
        T0(1)= l(T(j-1));
        T0(M)= r(T(j-1));
    end
    W = A*Wij(:,i-1)-Wij(:,i-2)+sigma^2*T0;
end


%%Solución real
ur = @(x,t) exp(-x-2*t);
u = zeros(M,N);
for i=1:M
    for j=1:N
        u(i,j) = ur(X(i),T(j));
    end
end
surf(X,T,u)
