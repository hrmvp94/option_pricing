S = 20;
E = 10;
r = 0.1;
sigma = 0.4;
T = 0.25;
M = 40;
Xmin = -10;
Xmax = 10;
dx = 1/10;
Nminus = Xmin/dx;
Nplus = Xmax/dx;
dt = 0.5*sigma^2*T/M;
values = zeros(Nplus - Nminus + 1, 1);
exact = option(log(S/E), 0.5*sigma^2*T);
values = implict_fd_1(values, dx, dt, M, Nminus, Nplus);

function u_val = option(x, r)
k = r/(0.5*0.4^2);
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

function v_l = implict_fd_1(values, dx, dt, M, Nminus, Nplus)
a = dt/(dx^2);
b = zeros(Nplus - Nminus + 1, 1);
for n = 1:1:Nplus - Nminus + 1
    values(n) = European_0(n*dx);
end
y = zeros(Nplus - Nminus + 1, 1);
y = lu_find_y(y, a, Nminus, Nplus);
for m = 1:1:M
    tau = m*dt;
    for n = 2:1:Nplus - Nminus
        b(n) = values(n);
    end
    values(1) = option(Nminus*dx, tau);
    values(Nplus - Nminus + 1) = option(Nplus*dx, tau);
    b(2) = b(2) + a*values(1);
    b(Nplus - Nminus) = b(Nplus - Nminus) + a*values(Nplus - Nminus + 1);
            
    values = lu_solver(values, b, y, a, Nminus, Nplus);
end
v_l = values;
end

function y_r_f = lu_find_y(y, a, Nminus, Nplus)
asq = a^2;
y(2) = 1 + 2*a;
for n = 3:1:Nplus - Nminus
    y(n) = 1 + 2*a - asq/y(n - 1);
end
y_r_f = y;
end

function u_r_l = lu_solver(u, b, y, a, Nminus, Nplus)
q = zeros(Nplus - Nminus + 1, 1);
q(2) = b(2);
for n = 3:1:Nplus - Nminus
    q(n) = b(n) + a*q(n - 1)/y(n - 1);
end
u(Nplus - Nminus) = q(Nplus - Nminus)/y(Nplus - Nminus);
for n = Nplus - Nminus - 1:-1:2
    u(n) = (q(n) + a*u(n + 1))/y(n);
end
u_r_l = u;
end


            

