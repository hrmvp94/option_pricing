S = 10;
E = 10;
r = 0.1;
sigma = 0.4;
T = 0.25;
M = 800;
Xmin = -10;
Xmax = 10;
dx = 1/100;
Nminus = Xmin/dx;
Nplus = Xmax/dx;
dt = 0.5*sigma^2*T/M;
alpha = dt/(dx^2);
L = eye(Nplus - Nminus - 1);
U = zeros(Nplus - Nminus - 1);
y_n = 1 + 2*alpha;
for n = 1:1:Nplus - Nminus - 1
    U(n, n) = y_n;
    if n < Nplus - Nminus - 1
        U(n, n + 1) = -alpha;
        L(n + 1, n) = -alpha/y_n;
    end
    y_n = (1 + 2*alpha) - alpha^2/y_n;
end
u_0 = zeros(Nplus - Nminus - 1, 1);
for n = 1:1:Nplus - Nminus - 1
    u_0(n) = European_0((Nminus + n)*dx);
end
u_0(1) = u_0(1) + alpha*option(Nminus*dx, dt);
u_0(Nplus - Nminus - 1) = u_0(Nplus - Nminus - 1) + alpha*option(Nplus*dx, dt);
q = L\u_0;
u_0 = U\q;
for m = 2:1:M
    u_0(1) = u_0(1) + alpha*option(Nminus*dx, m*dt);
    u_0(Nplus - Nminus - 1) = u_0(Nplus - Nminus - 1) + alpha*option(Nplus*dx, m*dt);
    q = L\u_0;
    u_0 = U\q;
end
exact = zeros(Nplus - Nminus - 1, 1);
for n = 1:1:Nplus - Nminus - 1
    exact(n) = option((Nminus + n)*dx, dt*M);
end
diff = norm(exact - u_0)
k = 0.1/(0.5*0.4^2);
true = zeros(Nplus - Nminus - 1, 1);
for m = 1:1:Nplus - Nminus - 1
    true(m) = exp(-0.5*(k - 1)*(Nminus + m)*dx - 0.25*(k + 1)^2*dt*M)*u_0(m)*E;
end

function u_val = option(x, r)
k = 0.1/(0.5*0.4^2);
d1 = x/sqrt(2*r) + 0.5*(k + 1)*sqrt(2*r);
I1 = exp(0.5*(k + 1)*x + 0.25*r*(k + 1)^2)*normcdf(d1);
d2 = x/sqrt(2*r) + 0.5*(k - 1)*sqrt(2*r);
I2 = exp(0.5*(k - 1)*x + 0.25*r*(k - 1)^2)*normcdf(d2);
u_val = I1 - I2;
end

function u_0 = European_0(x)
k = 0.1/(0.5*0.4^2);
payoff = exp(0.5*(k + 1)*x) - exp(0.5*(k - 1)*x);
if payoff >= 0
    u_0 = payoff;
else
    u_0 = 0;
end
end
