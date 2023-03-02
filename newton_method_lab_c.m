%%
psmax = 75000;
psmin = 100000;
ks = 0.045;
pfmax = 300000;
kf = 0.08;
p0 = 10000;

ps = @(t) psmax.*exp(-ks.*t) + psmin;
pf = @(t) pfmax/(1+((pfmax/p0)-1).*exp(-kf.*t));

%befolkningsökningen kan defineras som Pf - Ps
% funktion ska lösas för när förort är 20% större än stan (pf - 1,2ps = 0)

f = @(t) pfmax./(1+((pfmax./(p0))-1)*exp(-kf.*t)) - 1.2*(psmax*exp(-ks.*t)+psmin);
fp = @(t) 4050*exp(0.045.*t) + ((696000*exp(-0.08.*t))./((29*exp(-0.08.*t)+1).^2));

tol = 1e-8;
delta_t = 1;
t = 20;

guess = [];

while abs(delta_t) > tol
    tnext = t - f(t)/fp(t);
    delta_t = abs(tnext-t);
    t = tnext;
    guess = [guess, t];
end

diffT1 = abs(guess(2:end)-guess(1:end-1));
diffT2 = abs(diffT1(2:end)./(diffT1(1:end-1).^2));

kon = log2(abs(diffT2(2)-diffT2(1))/(abs(diffT2(3)-diffT2(2))).^2);

%%
% -----------------------------------2--------------------------------------

load("STOCKHOLMSDATA.mat");

c3 = [];
c1 = repelem(1, 70080);
c1 = c1';

c2 = 1:1:70080;
c2 = c2';

for n = 1:70080
    a = sin(((2*pi)/(365*24))*(n-(117*24)));
    c3 = [c3, a];
end

c3 = c3';

t = 0:1:70080;
t = t';

Ainterp = [c1 c2 c3];

A = Ainterp;
y = Tm;

cinterp = Ainterp\Tm;
x = 1:1:70080;
T = cinterp(1) + cinterp(2).*x + cinterp(3).*sin(((2*pi)/(365*24))*(x-(117*24)));
t240 = cinterp(1) + cinterp(2).*x + cinterp(3).*sin(((2*pi)/(365*24))*(x-240));

plot(x,Tm);

hold on

plot(x, t240, LineWidth=2);
plot(x,T, LineWidth=2)

res = Ainterp*cinterp - y;
resnorm = norm(res).^2;
disp(resnorm)

for i = 1:4368
    x = 0:0.1:100;
    t = cinterp(1) + cinterp(2).*x + cinterp(3).*sin(((2*pi)/(365*24))*(x-i));
    plot(x,t);
end

Tl = min(T);
Th = max(T);

yline(Tl, "g", 'DisplayName','Min');
yline(Th, "m", 'DisplayName','Max');

hold off

xlabel("Tid")
ylabel("Temperatur i C")

% -------------------------------------3------------------------------------
%%
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

disp("The integral is: " + TSim)
%%
% -------------------------------------4------------------------------------
%
t0 = 0;
tspan = [t0, 23];
y0 = [0 0];
[TOUT, YOUT]= ode45(@odesyst, tspan, y0);
y = YOUT(1:end,1);
ymax = max(y);

if ymax > 100
    disp("VARNING: DU KOMMER SLÅ I BACKEN")
else
    disp("DU KOMMER INTE SLÅ I BACKEN")
end
