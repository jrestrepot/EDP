M = 20;
N = 20;
a = 0;
b = 1;
at = 0;
bt = 1;
dd = @(x,t)4*exp(-x-2*t); %Segunda derivdada de u
f = @(x)exp(-x); %Condicion inicial
g= @(x)-2*exp(-x); %Condicion inicial
l = @(t)exp(-2*t); %Condicion de frontera
r = @(t)exp(-1-2*t); %Condicion de frontera
m =M+1; n=N+1; mn = m*n;
h = (b-a)/M; 
k =(bt-at)/4; 
x = a+(0:M)*h;
t = at+(0:N)*k;
A = zeros(mn,mn); 
b= zeros(mn,1);

for i=1:m
    j =1;
    A(i+(j-1)*m, i+(j-1)*m) = 1;
    b(i+(j-1)*m) = f(x(i));
    j=n;
    A(i+(j-1)*m, i+(j-1)*m)=1;
    b(i+(j-1)*m) = g(x(i));
end
for j=2:n-1
    i =1;
    A(i+(j-1)*m, i+(j-1)*m) = 1;
    b(i+(j-1)*m) = l(t(j));
    i = m;
    A(i+(j-1)*m, i+(j-1)*m) = 1;
    b(i+(j-1)*m) = r(t(j));
end

for i=2:m-1
    for j=2:n-1
     A(i+(j-1)*m, i-1+(j-1)*m) = 1/(h^2);
     A(i+(j-1)*m, i+1+(j-1)*m) = 1/(h^2);
     A(i+(j-1)*m, i+(j-1)*m) = -3/h^2-2/k^2;
     A(i+(j-1)*m, i+(j-2)*m) = 1/(k^2);
     A(i+(j-1)*m, i+j*m) = 1/(k^2);
     b(i+(j-1)*m)=dd(x(i),t(j));
    end
end
v=A\b
w=reshape(v(1:mn),m,n);
surf(x,t,w)

%%Solución real
ur = @(x,t) exp(-x-2*t);
u = zeros(m,n);
for i=1:m
    for j=1:n
        u(i,j) = ur(x(i),t(j));
    end
end
surf(x,t,u)
