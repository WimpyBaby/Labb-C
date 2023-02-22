% Q = @(t) (9+5*(cos(0.2*t)).^2);
% c = @(t) (5*exp(-0.5*t)+2*exp(0.15*t));
% f = @(t) Q(t).*c(t);
% 
% Shvec = [];
% 
% for n = [10 20 40]
%     h = 6/n;
%     t = 3:h:9;
%     Sh = ((h/3)*(sum(f(t(2:2:end-1)))*4+(sum(f(t(3:2:end-2)))*2+f(t(1))+f(t(end)))));
%     Shvec = [Shvec, Sh];
% end
% 
% nog = log2(abs(Shvec(2)-Shvec(1))/abs(Shvec(3)-Shvec(2)));

% psmax = 75000;
% psmin = 100000;
% ks = 0.045;
% pfmax = 300000;
% kf = 0.08;
% p0 = 10000;
% 
% ps = @(t) psmax.*exp(-ks.*t) + psmin;
% pf = @(t) pfmax/(1+((pfmax/p0)-1).*exp(-kf.*t));
% 
% f = @(t) pfmax./(1+((pfmax./(p0))-1)*exp(-kf.*t)) - 1.2*(psmax*exp(-ks.*t)+psmin);
% fp = @(t) 4050*exp(0.045.*t) + ((696000*exp(-0.08.*t))./((29*exp(-0.08.*t)+1).^2));
% 
% t = 35;
% tol = 1e-8;
% diff = 1
% 
% tvec = [ t ];
% 
% while abs(diff) > tol
%     tny = t - f(t)/fp(t);
%     diff = abs(tny - t);
%     t = tny;
%     tvec = [tvec, t];
% end
% 
% diffT1 = abs(tvec(2:end)-tvec(1:end-1));
% diffT2 = abs(diffT1(2:end)./diffT1(1:end-1).^2);
% diffT2 = diffT2';

h = 0.25;

x = [0 0.125 0.25 0.375 0.5 0.625 0.75 0.875 1];
b = [2.6 2.08 1.72 1.45 1.26 1.13 1.04 0.97 0.92];
f = b.*x;

T2h = (h/2).*(f(1) + 2.*f(3) + 2.*f(5) + 2.*f(7) + f(9));
