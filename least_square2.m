load("STOCKHOLMSDATA.mat");

w = (2*pi)/(365*24);
t = (1:70080)';
col1 = ones(70080, 1);
col2 = t;
col240 = sin(w*(t-240));

A240 = [col1, col2, col240];
c240 = A240\Tm;
T240 = c240(1) + c240(2)*t + c240(3)*sin(w*(t-240));

plot(t, Tm);

hold on

plot(t, T240)
findts = [zeros(4368, 1)]; 

resnorm = [];

for ts = 1:4368

    col3 = sin(((2*pi)/(365*24))*(t-ts));
    A = [col1, col2, col3];
    c = A\Tm;
    T = c(1) + c(2)*t + c(3)*sin(w*(t-ts));
    res = Tm - A*c;
    resnorm = [resnorm, norm(res).^2];
end

minres = min(resnorm);



