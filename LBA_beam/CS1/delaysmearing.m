function att = delaysmearing(l, df, u0)

dr = 0.1;
da = 0.01;
c = 2.99792e8;
const = 2 * pi * df / c;
r = 0:dr:6 * u0;
a = 0:da:const * max(r) * max(l);
psi = 0:0.01:2*pi;

disp('doing integral over psi');
tic
fa = zeros(length(a), 1);
for idx = 2:min(length(a), 1000)
    fa(idx) = sum(sin(a(idx) * cos(psi)) ./ (a(idx) * cos(psi)) ) * da;
end
fa(1) = 2 * pi;
if (length(a) > 1000)
    fa(1001:length(a)) = 2 * pi ./ a(1001:end);
end
toc

disp('doing integral over r');
tic
att = zeros(length(l), 1);
for idx = 1:length(l);
    att(idx) = sum(interp1(a, fa, const * r * l(idx)) .* exp(-0.5 * (r / u0).^2));
end
att = att / (u0 * pi * sqrt(2 * pi)) * dr;
toc
