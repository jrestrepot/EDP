%Código basado en el esquema matricial
c = 2;
m = 20; %#intervalos de x
f = @(x)exp(-x); &%cond inicial
g = @(x)-2*exp(-x); %cond inicial
l = @(t)exp(-2*t); %cond frontera
r = @(t) exp(-1-2*t); %cond frontera
h = 1/m; %definición de h
k = h/(c*2);  %garantiza estabilidad
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
   Fx(i,1) = f(X(i)); %sols de fx
   Gx(i,1) = g(X(i));  %sols de gx
end

W0 = Fx;
W1 = (1/2)*A*Fx +k*Gx+(sigma^2)/2*T0; %Primer paso temporal
Wij = ones(M,N);
Wij(:,1) = W0; %Se define la primera columna de la matriz como W0
Wij(:,2) = W1; %Se define la segunda columna de la matriz como W1

for i=3:N
    for j=2:M
        T0(1)= l(T(j-1));
        T0(M)= r(T(j-1));
        Wij(:,i) = A*Wij(:,i-1)-Wij(:,i-2)+sigma^2*T0; %matriz de soluciones, se llena por columnas
    end   
end

%Si se quieren ver las soluciones imprima Wij

%%Solución real
ur = @(x,t) exp(-x-2*t);
u = zeros(M,N);
for i=1:M
    for j=1:N
        u(i,j) = ur(X(i),T(j));
    end
end
surf(X,T,u)
