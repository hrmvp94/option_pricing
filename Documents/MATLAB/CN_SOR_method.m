S = 10;
E = 10;
r = 0.1;
sigma = 0.4;
T = 0.5;
M = 100;
Xmin = -10;
Xmax = 10;
dx = 1/50;
Nminus = Xmin/dx;
Nplus = Xmax/dx;
dt = 0.5*sigma^2*T/M;
alpha = dt/(dx^2);
eps = 10^(-9);
omega = 1;
domega = 0.05;
oldloops = 10000;
u_0 = zeros(Nplus - Nminus - 1, 1);
for n = 1:1:Nplus - Nminus - 1
    u_0(n) = European_0((Nminus + n)*dx);
end
b_n = zeros(Nplus - Nminus - 1, 1);
u_m_inf = European_0(Nminus*dx);
u_p_inf = European_0(Nplus*dx);
for m = 1:1:M
    tau = m*dt;
    error = 1;
    loops = 0;
    for n = 1:1:Nplus - Nminus - 1
        if n == 1
           b_n(n) = (1 - alpha)*u_0(n) + 0.5*alpha*(u_m_inf + u_0(n + 1));
        elseif n == Nplus - Nminus - 1
           b_n(n) = (1 - alpha)*u_0(n) + 0.5*alpha*(u_0(n - 1) + u_p_inf);
        else
           b_n(n) = (1 - alpha)*u_0(n) + 0.5*alpha*(u_0(n - 1) + u_0(n + 1));
        end
    end
    u_m_inf = option(Nminus*dx, tau);
    u_p_inf = option(Nplus*dx, tau);
    while error > eps
        error = 0;
        for n = 1:1:Nplus - Nminus - 1
           if n == 1
                y = (b_n(n) + 0.5*alpha*(u_m_inf + u_0(n + 1)))/(1 + alpha);
            elseif n == Nplus - Nminus - 1
                y = (b_n(n) + 0.5*alpha*(u_0(n - 1) + u_p_inf))/(1 + alpha);
            else
                y = (b_n(n) + 0.5*alpha*(u_0(n - 1) + u_0(n + 1)))/(1 + alpha);
           end
            y = u_0(n) + omega*(y - u_0(n));
            error = error + (u_0(n) - y)^2;
            u_0(n) = y;
        end
        loops = loops + 1;
    end
    if loops > oldloops
       domega = -domega;
    end
    omega = omega + domega;
    oldloops = loops;
end
exact = zeros(Nplus - Nminus - 1, 1);
for n = 1:1:Nplus - Nminus - 1
    exact(n) = option((Nminus + n)*dx, dt*M);
end
%diff = norm(exact - u_0)
k = 0.1/(0.5*0.4^2);
true = zeros(Nplus - Nminus - 1, 1);
for m = 1:1:Nplus - Nminus - 1
    true(m) = exp(-0.5*(k - 1)*(Nminus + m)*dx - 0.25*(k + 1)^2*dt*M)*u_0(m)*E;
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
