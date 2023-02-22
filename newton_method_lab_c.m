clear all;

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

% diffT1 = abs(guess(2:end)-guess(1:end-1));
% diffT2 = abs(diffT1(2:end)./(diffT1(1:end-1).^2));

kon = log2(abs(guess(2)-guess(1))/(abs(guess(3)-guess(2))).^2);

% --------------------------------------------------------------------------

% load("STOCKHOLMSDATA.mat");
% 
% c3 = [];
% c1 = repelem(1, 70080);
% c1 = c1';
% 
% c2 = 0:1:70079;
% c2 = c2';
% %c3 = sin(2.*pi./(365*24).*((0:1:70080);
% 
% for n = 0:70079
%     a = sin(((2*pi)/(365*24))*(n-(117*24)));
%     c3 = [c3, a];
% end
% 
% c3 = c3';
% 
% t = 0:1:70080;
% t = t';
% 
% Ainterp = [c1 c2 c3];
% 
% cinterp = Ainterp\Tm;
% x = 0:1:70079;
% % T = cinterp(1) + cinterp(2).*x + cinterp(3).*sin(((2*pi)/(365*24))*(x-240));
% T = @(x) cinterp(1) + cinterp(2).*x + cinterp(3).*sin(((2*pi)/(365*24))*(x-(117*24)));
% % T = T';
% res = Ainterp * cinterp - Tm;
% resnorm = norm(res).^2;
% disp(resnorm)


% % x = 0:1:70080;
% % i = lsqr(Ainterp, Tm);
% % T = i(1) + i(2) + i(3).*sin(((2*pi)/(365*24))*(x-240));

% plot(x, T(x))

% xlabel("Tid")
% ylabel("Temperatur i C")

% %T = @(t) c(1) + c(2).*t + c(3).*sin(2.*pi.*(t-240));
% 
% %Tmc = y

% --------------------------------------------------------------------------

% Q = @(t) (9+5*(cos(0.2*t)).^2);
% c = @(t) (5*exp(-0.5*t)+2*exp(0.15*t));
% f = @(t) Q(t).*c(t);
% 
% x0 = 3;
% x = 9;

% n = [10 20 40];
% n = 10;

% h = (x-x0)/n;
% t = x0:h:x;
% t2 = x0:0.3:x;
% t4 = x0:0.15:x;

% intvalues = [];
% 
% for i = n
%     h = (x-x0)/i;
%     t = x0:h:x;
%     TSim = ((h/3)*(sum(f(t(2:2:end-1)))*4+(sum(f(t(3:2:end-2)))*2+f(t(1))+f(t(end)))));
%     intvalues = [intvalues, TSim];
% end

% Tsim = (h/3).*(f(t(1))+(sum(f(t(2:2:end-1)))*4 +(sum(f(t(3:2:end-2)))*2 + f(t(end)))));
% T2Sim = (h/6).*(f(t2(1))+(sum(f(t2(2:2:end-1)))*4 +(sum(f(t2(3:2:end-2)))*2 + f(t2(end)))));
% T4Sim = (h/12).*(f(t4(1))+(sum(f(t4(2:2:end-1)))*4 +(sum(f(t4(3:2:end-2)))*2 + f(t4(end)))));
% 
% intvalues = [Tsim, T2Sim, T4Sim];
% 
% nog = log2(abs(intvalues(2)-intvalues(1))/abs(intvalues(3)-intvalues(2)));
