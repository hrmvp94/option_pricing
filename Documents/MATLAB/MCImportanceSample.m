randn('seed', 0);
S0 = 50;
K = 52;
r = 0.1;
T = 2/12;
sigma = 0.4;
Sb = 30;
NSteps = 60;
NRepl = 10000;
bp = 200;
[P, CI, NCrossed] = DOPutMCCondIS(S0, K, r, T, sigma, Sb, NSteps, NRepl, bp)
function [Pdo, CI, NCrossed] = DOPutMCCondIS(S0, K, r, T, sigma, Sb, NSteps, NRepl, bp)
dt = T/NSteps;
nudt = (r - 0.5*sigma^2)*dt;
b = bp*nudt;
sidt = sigma*sqrt(dt);
[Call, Put] = BS(S0, K, r, T, sigma);
NCrossed = 0;
Payoff = zeros(NRepl, 1);
Times = zeros(NRepl, 1);
StockVals = zeros(NRepl, 1);
ISRatio = zeros(NRepl, 1);
for i = 1:NRepl
    vetZ = nudt - b + sidt*randn(1, NSteps);
    LogPath = cumsum([log(S0), vetZ]);
    Path = exp(LogPath);
    jcrossed = min(find(Path <= Sb));
    if not(isempty(jcrossed))
        NCrossed = NCrossed + 1;
        TBreach = jcrossed - 1;
        Times(NCrossed) = TBreach*dt;
        StockVals(NCrossed) = Path(jcrossed);
        ISRatio(NCrossed) = exp(TBreach*b^2/2/sigma^2/dt + ...
            b/sigma^2/dt*sum(vetZ(1:TBreach)) - ...
            TBreach*b/sigma^2*(r - sigma^2/2));
    end
end
if (NCrossed > 0)
    [Caux, Paux] = BS(StockVals(1:NCrossed), K, r, ...
        T - Times(1:NCrossed), sigma);
    Payoff(1:NCrossed) = exp(-r*Times(1:NCrossed)).* Paux ...
        .* ISRatio(1:NCrossed);
end
[Pdo, aux, CI] = normfit(Put - Payoff);
end

function [Call, Put] = BS(S, E, r, T, vol)
sz = size(S, 1);
Call = zeros(sz, 1);
Put = zeros(sz, 1);
for i = 1:1:sz
    d1 = (log(S(i)/E) + (r + 0.5*vol^2)*T(i))/(vol*sqrt(T(i)));
    d2 = d1 - vol*sqrt(T(i));
    Call(i) = S(i)*normcdf(d1) - E*exp(-r*T(i))*normcdf(d2);
    Put(i) = E*exp(-r*T(i))*normcdf(-d2) - S(i)*normcdf(-d1);
end
end