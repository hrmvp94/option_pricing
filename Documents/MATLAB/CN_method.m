S = 20;
E = 10;
r = 0.1;
sigma = 0.4;
T = 0.25;
M = 40;
Xmin = -10;
Xmax = 10;
dx = 1/40;
Nminus = Xmin/dx;
Nplus = Xmax/dx;
dt = 0.5*sigma^2*T/M;
alpha = dt/(dx^2);
A = gallery('tridiag',Nplus - Nminus - 1, -0.5*alpha,1 + alpha,-0.5*alpha);
u_0 = zeros(Nplus - Nminus - 1, 1);
for n = 1:1:Nplus - Nminus - 1
    u_0(n) = European_0((Nminus + n)*dx);
end
z_0 = zeros(Nplus - Nminus - 1, 1);
for n = 1:1:Nplus - Nminus - 1
    if n == 1
        z_0(n) = (1 - alpha)*u_0(n) + 0.5*alpha*(option(Nminus*dx, 0) + u_0(n + 1));
    elseif n == Nplus - Nminus - 1
        z_0(n) = (1 - alpha)*u_0(n) + 0.5*alpha*(u_0(n - 1) + option(Nplus*dx, 0));
    else
        z_0(n) = (1 - alpha)*u_0(n) + 0.5*alpha*(u_0(n - 1) + u_0(n + 1));
    end
end
z_0(1) = z_0(1) + 0.5*alpha*option(Nminus*dx, dt);
z_0(Nplus - Nminus - 1) = z_0(Nplus - Nminus - 1) + 0.5*alpha*option(Nplus*dx, dt);
u_0 = A\z_0;
for m = 1:1:M - 1
    for n = 1:1:Nplus - Nminus - 1
        if n == 1
            z_0(n) = (1 - alpha)*u_0(n) + 0.5*alpha*(option(Nminus*dx, m*dt) + u_0(n + 1));
        elseif n == Nplus - Nminus - 1
            z_0(n) = (1 - alpha)*u_0(n) + 0.5*alpha*(u_0(n - 1) + option(Nplus*dx, m*dt));
        else
            z_0(n) = (1 - alpha)*u_0(n) + 0.5*alpha*(u_0(n - 1) + u_0(n + 1));
        end
    end
    z_0(1) = z_0(1) + 0.5*alpha*option(Nminus*dx, (m + 1)*dt);
    z_0(Nplus - Nminus - 1) = z_0(Nplus - Nminus - 1) + 0.5*alpha*option(Nplus*dx, (m + 1)*dt);
    u_0 = A\z_0;
end
exact = zeros(Nplus - Nminus - 1, 1);
for n = 1:1:Nplus - Nminus - 1
    exact(n) = option((Nminus + n)*dx, dt*M);
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