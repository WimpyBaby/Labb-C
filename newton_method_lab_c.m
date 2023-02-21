clear all;

% psmax = 75000;
% psmin = 100000;
% ks = 0.045;
% pfmax = 300000;
% kf = 0.08;
% p0 = 10000;
% 
% ps = @(t) psmax.*exp(-ks.*t) + psmin;
% pf = @(t) pfmax/(1+((pfmax/p0)-1).*exp(-kf.*t));

%befolkningsökningen kan defineras som Pf - Ps
% funktion ska lösas för när förort är 20% större än stan (pf - 1,2ps = 0)

% f = @(t) pfmax./(1+((pfmax./(p0))-1)*exp(-kf.*t)) - 1.2*(psmax*exp(-ks.*t)+psmin);
% fp = @(t) 4050*exp(0.045.*t) + ((696000*exp(-0.08.*t))./((29*exp(-0.08.*t)+1).^2));
% 
% tol = 1e-8;
% delta_t = 1;
% t = 20;
% 
% guess = [];
% 
% while delta_t > tol
%     tnext = t - f(t)/fp(t);
%     delta_t = abs(tnext-t);
%     t = tnext;
%     guess = [guess, t];
% end
% 
% Hn = abs(guess(2:end)-guess(1:end-1));
% konv = 

% --------------------------------------------------------------------------

load("STOCKHOLMSDATA.mat");

c3 = [];
c1 = repelem(1, 70080);
c1 = c1';

c2 = 0:1:70079;
c2 = c2';
%c3 = sin(2.*pi./(365*24).*((0:1:70080);

for n = 0:70079
    a = sin(((2*pi)/(365*24))*(n-240));
    c3 = [c3, a];
end

c3 = c3';

t = 0:1:70080;
t = t';

Ainterp = [c1 c2 c3];

cinterp = Ainterp\Tm;
x = 0:1:70079;
T = cinterp(1) + cinterp(2).*x + cinterp(3).*sin(((2*pi)/(365*24))*(x-240));
T = T';
residual = abs(norm((Tm-T).^2));

% x = 0:1:70080;
% i = lsqr(Ainterp, Tm);
% T = i(1) + i(2) + i(3).*sin(((2*pi)/(365*24))*(x-240));

plot(x, T)

% xlabel("Tid i timmar")
% ylabel("Temperatur i C")

%T = @(t) c(1) + c(2).*t + c(3).*sin(2.*pi.*(t-240));

%Tmc = y




