Q = @(t) (9+5*(cos(0.4*t)).^2);
c = @(t) (5*exp(-0.5*t)+2*exp(0.15*t));
f = @(t) Q(t).*c(t);

x0 = 3;
x = 9;
n = [10 20 40];


intvalues = [];

for i = n
    h = (x-x0)/i;
    t = x0:h:x;
    TSim = (h/3)*(4*(sum(f(t(2:2:end-1)))) + 2*(sum(f(t(3:2:end-2)))) + f(t(1)) + f(t(end)));
    intvalues = [intvalues, TSim];
end

nog = log2(abs(intvalues(2)-intvalues(1))/abs(intvalues(3)-intvalues(2)));

disp("Massan är : " + TSim)