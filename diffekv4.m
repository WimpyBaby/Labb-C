y0 = [0;0];
tspan = [0 60];

[t,y] = ode45(@odefun, tspan, y0);

ymax = max(y);

if ymax < 100
    disp("You are well")
else
    disp("not doing good")
end

plot(t, y(:,1))

function f=odefun(t,y)
g = 9.81;
cd = 0.25;
k = 40;
lambda = 8;
L = 30;
m = 65;

u = y(1);
v = y(2);

if u <= 30
    dv = g - sign(v)*((cd/m)*v)^2;
    f = [v; dv];
else
    dv = g - sign(v)*((cd/m)*v)^2 - (k/m)*(u-L)-(lambda/m)*v;
    f = [v; dv];
end
end