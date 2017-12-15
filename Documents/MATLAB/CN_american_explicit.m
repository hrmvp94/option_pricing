S = 10;
E = 10;
r = 0.1;
sigma = 0.4;
T = 0.5;
M = 1024;
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
g_n = zeros(Nplus - Nminus + 1, 1);
for m = 1:1:M
    q(1) = option(Nminus*dx, m*dt);
    q(Nplus - Nminus + 1) = option(Nplus*dx, m*dt);
     for n = 1:1:Nplus - Nminus + 1
        g_n(n) = option((Nminus + n - 1)*dx, m*dt);
     end
    for n = 2:1:Nplus - Nminus
        q(n) = alpha*u_0(n - 1) + (1 - 2*alpha)*u_0(n) + alpha*u_0(n + 1);
        u_0(n) = max(q(n), g_n(n));
    end
    u_0(1) = q(1);
    u_0(Nplus - Nminus + 1) = q(Nplus - Nminus + 1);
end
k = 0.1/(0.5*0.4^2);
true = zeros(Nplus - Nminus + 1, 1);
%ime = zeros(Nplus - Nminus - 1, 1);
for m = 1:1:Nplus - Nminus + 1
    true(m) = exp(-0.5*(k - 1)*(Nminus + m - 1)*dx - 0.25*(k + 1)^2*dt*M)*u_0(m)*E;
    %ime(m) = exp(-0.5*(k - 1)*(Nminus + m)*dx - 0.25*(k + 1)^2*dt*M)*g_n(m)*E;
end


function u_val = option(x, r)
k = 0.1/(0.5*0.4^2);
d1 = x/sqrt(2*r) + 0.5*(k + 1)*sqrt(2*r);
I1 = exp(0.5*(k + 1)*x + 0.25*r*(k + 1)^2)*normcdf(-d1);
d2 = x/sqrt(2*r) + 0.5*(k - 1)*sqrt(2*r);
I2 = exp(0.5*(k - 1)*x + 0.25*r*(k - 1)^2)*normcdf(-d2);
u_val = I2 - I1;
end

function u_0 = European_0(x)
k = 0.1/(0.5*0.4^2);
payoff = exp(0.5*(k - 1)*x) - exp(0.5*(k + 1)*x);
if payoff >= 0
    u_0 = payoff;
else
    u_0 = 0;
end
end