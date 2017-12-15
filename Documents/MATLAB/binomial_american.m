S = 10;
E = 10;
r = 0.1;
sigma = 0.4;
T = 0.25;
M = 128;
dt = T/M;
discount = exp(-r*dt);
A = 0.5*(discount + exp((r + sigma^2)*dt));
u = A + sqrt(A^2 - 1);
d = A - sqrt(A^2 - 1);
p = (exp(r*dt) - d)/(u - d);
S_0 = zeros(M + 1, M + 1);
for i = 1:1:M + 1
    for j = i:1:M + 1
        S_0(i, j) = S*u^(j - i)*d^(i - 1);
    end
end
v = zeros(M + 1, M + 1);
for i = 1:1:M + 1
    v(i, M + 1) = payoff(S_0(i, M + 1), E);
end
for j = M:-1:1
    for i = 1:1:j
        hold = discount*(p*v(i, j + 1) + (1 - p)*v(i + 1, j + 1));
        v(i, j) = max(hold, payoff(S_0(i, j), E));
    end
end
diff = norm(0.6913 - v(1, 1))

function p = payoff(S, K)
if S < K
    p = K - S;
else
    p = 0;
end
end