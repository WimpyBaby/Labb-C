clear all;

psmax = 75000;
psmin = 100000;
ks = 0.045;
pfmax = 300000;
kf = 0.08;
p0 = 10000;

ps = @(t) psmax*exp(-ks.*t)+psmin;
pf = @(t) pfmax./(1+((pfmax./(p0))-1)*exp(-kf.*t));

f = @(t) pfmax./(1+((pfmax./(p0))-1)*exp(-kf.*t)) - 1.2*(psmax*exp(-ks.*t)+psmin);
fprim = @(t) 4050*exp(0.045.*t) + ((696000*exp(-0.08.*t))./((29*exp(-0.08.*t)+1).^2));

t = 0:1:100;
plot(t, f(t))

t = 35;
tol = 1e-8;
diff = 1;

tvec = [ t ];

while abs(diff)>tol
    tny = t - f(t)/fprim(t);
    diff = abs(tny-t);
    t = tny;
    tvec = [tvec, t];
end

disp(["skiten blir: " num2str(t)])