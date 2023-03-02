psmax = 75000;
psmin = 100000;
ks = 0.045;
pfmax = 300000;
kf = 0.08;
p0 = 10000;

f = @(t) pfmax./(1+((pfmax./(p0))-1)*exp(-kf.*t)) - 1.2*(psmax*exp(-ks.*t)+psmin);
fp = @(t) 4050*exp(0.045.*t) + ((696000*exp(-0.08.*t))./((29*exp(-0.08.*t)+1).^2));

t = 20;
guess = [];
tol = 1e-4;
iter = 100;
diff = zeros(iter, 1);

for i = 1:iter
    tnext = t - f(t)/fp(t);
    diff(i) = abs(t-tnext);
    t = tnext;
    if diff < tol
        break
    end
    n = length(diff);
    iter = 1:n;
    guess = [guess, t];
    loglog(iter, diff, '-o')
end

%diff = 1;

% while diff > tol
%     tnext = t - f(t)/fp(t);
%     diff = abs(t-tnext);
%     t = tnext;
%     guess = [guess, t];
% end

% xdiff = abs(guess(2:end)-guess(1:end-1));
% k = abs(xdiff(2:end)) ./ (xdiff(1:end-1)).^2;

