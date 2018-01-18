%S = 0;
E = 10;
r = 0.1;
sigma = 0.4;
T = 0.5;
N = 400;
Xmin = 0;
Xmax = 50;
dx = 1/16;
Nminus = Xmin/dx;
Nplus = Xmax/dx;
%dx = log(S/E)/1280;
dt = T/N;
alpha = dt/(dx^2);
eps = 10^(-9);

u_0 = zeros(Nplus - Nminus - 1, 1);
for n = 1:1:Nplus - Nminus - 1
    u_0(n) = BS(E, (Nminus + n)*dx, r, sigma,0, 0, T, 'put');
    %u_0(n) = European_0((Nminus + n)*dx);
end
g_n = zeros(Nplus - Nminus - 1, 1);
u_m_inf = BS(E, Nminus*dx, r, sigma,0, 0, T, 'put');
u_p_inf = BS(E, Nplus*dx, r, sigma,0, 0, T, 'put');
M = gallery('tridiag',Nplus - Nminus - 1,alpha,-2*alpha,alpha);
I = eye(Nplus - Nminus - 1);
for m = 1:1:N
    tau = m*dt;
    initial = u_0;
    error = 1;
    for n = 1:1:Nplus - Nminus - 1
        g_n(n) = g_payoff((Nminus + n)*dx, E);
    end
    u_m_inf = BS(E, Nminus*dx, r, sigma,0, tau, T, 'put');
    u_p_inf = BS(E, Nplus*dx, r, sigma,0, tau, T, 'put');
    P_prev = I;
    P_cur = I;
    for k = 1:1:Nplus - Nminus - 1
        if initial(k) < g_n(k)
            P_cur(k, k) = 1/eps;
        else
            P_cur(k, k) = 0;
        end
    end
    while (~isequal(P_prev, P_cur)) && (error > 10^(-8))
        tic
        %F = (I - 0.5*M)*initial - (I + 0.5*M)*u_0 + P_cur*(initial - g_n);
        F = (I - M)*initial - I*u_0 + P_cur*(initial - g_n);
        %F(1) = F(1) - (option(Nminus*dx, tau) - option(Nminus*dx, tau - dt))*0.5*alpha;
        %F(Nplus - Nminus - 1) = F(Nplus - Nminus - 1) - (option(Nplus*dx, tau) - option(Nplus*dx, tau - dt))*0.5*alpha;
        %F_p = I - 0.5*M + P_cur;
        F(1) = F(1) - u_m_inf*alpha;
        F(Nplus - Nminus - 1) = F(Nplus - Nminus - 1) - u_p_inf*alpha;
        F_p = I - M + P_cur;
        initial_prev = initial;
        initial = initial - F_p\F;
        toc
        P_prev = P_cur;
        for j = 1:1:Nplus - Nminus - 1
            if initial(j) < g_n(j)
                P_cur(j, j) = 1/eps;
            else
                P_cur(j, j) = 0;
            end
        end
        error = norm((initial - initial_prev)/max(1, norm(initial)), 1);
    end
    u_0 = initial;
end

%slope = (u_0(4895) - u_0(4894))/dx;
%intercept = u_0(4894) - (Xmin + 4894*dx)*slope;
%value = log(S/E)*slope + intercept;
% i_true = zeros(Nplus - Nminus - 1, 1);
% axis = zeros(Nplus - Nminus - 1, 1);
% for m = 1:1:Nplus - Nminus - 1
%     i_true(m) = exp(-0.5*(k - 1)*(Nminus + m)*dx - 0.25*(k + 1)^2*dt*N)*u_0(m)*E;
%     axis(m) = E*exp((Nminus + m)*dx);
% end

function b = BS(E,S, r, vol,d, t, T, option)
T1 = T - t;
d1 = (log(S/E) + (r - d + 0.5*vol^2)*T1)/(vol*sqrt(T1));
d2 = d1 - vol*sqrt(T1);
switch option
    case 'call'
        b = S*exp(-d*T1)*normcdf(d1) - E*exp(-r*T1)*normcdf(d2);
    case 'put'
        b = E*exp(-r*T1)*normcdf(-d2) - S*exp(-d*T1)*normcdf(-d1);
end
end

function g_p = g_payoff(S, E)
payoff = E - S;
if payoff >= 0
    g_p = payoff;
else
    g_p = 0;
end
end