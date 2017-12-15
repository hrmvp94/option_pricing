S = 20;
E = 10;
r = 0.1;
sigma = 0.4;
T = 0.25;
M = 32;
Xmin = -10;
Xmax = 10;
dx = 1/80;
Nminus = Xmin/dx;
Nplus = Xmax/dx;
dt = 0.5*sigma^2*T/M;
alpha = dt/(dx^2);
u_0 = zeros(Nplus - Nminus + 1, 1);
for n = 1:1:Nplus - Nminus + 1
    u_0(n) = European_0((Nminus + n - 1)*dx);
end
q = zeros(Nplus - Nminus + 1, 1);
for m = 1:1:M
    q(1) = option(Nminus*dx, m*dt);
    q(Nplus - Nminus + 1) = option(Nplus*dx, m*dt);
    for n = 2:1:Nplus - Nminus
        q(n) = alpha*u_0(n - 1) + (1 - 2*alpha)*u_0(n) + alpha*u_0(n + 1);
    end
    u_0 = q;
end
exact = zeros(Nplus - Nminus + 1, 1);
for n = 1:1:Nplus - Nminus + 1
    exact(n) = option((Nminus + n - 1)*dx, dt*M);
end
diff = norm(exact - u_0)


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