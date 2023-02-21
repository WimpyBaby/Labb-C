h1 = 0.25;
h2 = 0.125;

x = 0:0.125:1;

f = [2.6 2.08 1.72 1.45 1.26 1.13 1.04 0.97 0.92];

T2h = (h1/2).*(f(1).*x(1) + 2.*(f(3).*x(3) + f(5).*x(5) + f(7).*x(7)) + f(9)*x(9)); 

Th = (h2/2).*(f(1).*x(1) + 2.*(f(2).*x(2) + f(3).*x(3) + f(4).*x(4) + f(5).*x(5) + f(6).*x(6) + f(7).*x(7) + f(8).*x(8)) + f(9).*x(9));

Tstar = (4.*Th - T2h)./3;

disp(Tstar)