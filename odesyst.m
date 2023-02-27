function f=odesyst(t,y)
g = 9.81;
cd = 0.25;
k = 40;
lambda = 8;
L = 30;
m = 60;

u = y(1);
v = y(2);

if u<=30
    du = v;
    dv = g-sign(v)*(cd/m)*(v^2);
    f = [du; dv];
else u > 30;
    du = v;
    dv = g-sign(v)*(cd/m)*(v^2)-((k/m)*(u-L))-((lambda/m)*v);
    f = [du; dv];
end

end