randn('state', 0);
NRepl = 1000;
T = 1;
NSteps = 4;
WSamples = zeros(NRepl, 1 + NSteps);
for i = 1:NRepl
    WSamples(i,:) = WienerBridge(T, NSteps);
end
m = mean(WSamples(:, 2:(1 + NSteps)));
sd = sqrt(var(WSamples(:, 2:(1 + NSteps))));

GBMBridge(50, 0.5, 0.4, 1, 64, 20)

function SPaths = GBMBridge(S0, mu, sigma, T, NSteps, NRepl)
dt = T/NSteps;
nudt = (mu - 0.5*sigma^2)*dt;
SPaths = zeros(NRepl, NSteps + 1);
for k = 1:NRepl
    W = WienerBridge(T, NSteps);
    Increments = nudt + sigma*diff(W');
    LogPath = cumsum([log(S0), Increments]);
    SPaths(k, :) = exp(LogPath);
end
SPaths(:, 1) = S0;
end

function WSamples = WienerBridge(T, NSteps)
NBisections = log2(NSteps);
if round(NBisections) ~= NBisections
    fprintf('ERROR: NSteps must be a power of 2\n');
    return
end
WSamples = zeros(NSteps + 1, 1);
WSamples(1) = 0;
WSamples(NSteps + 1) = sqrt(T)*randn;
TJump = T;
IJump = NSteps;
for k = 1:NBisections
    left = 1;
    i = IJump/2 + 1;
    right = IJump + 1;
    for j = 1:2^(k - 1)
        a = 0.5*(WSamples(left) + WSamples(right));
        b = 0.5*sqrt(TJump);
        WSamples(i) = a + b*randn;
        right = right + IJump;
        left = left + IJump;
        i = i + IJump;
    end
    IJump = IJump/2;
    TJump = TJump/2;
end
end