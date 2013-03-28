% core
xpos = [0, sin(2 * pi * (0:4) / 5) * 2.4];
ypos = [0, cos(2 * pi * (0:4) / 5) * 2.4];

% rings
radius = [5.8 9.7 14.1 19.5 25.9 33.3 41.7];
nelem = [9 13 17 11 13 15 9];
for idx = 1:length(radius)
  disp(num2str(idx));
  xpos = [xpos, sin(2 * pi * (1:nelem(idx)) / nelem(idx)) * radius(idx)];
  ypos = [ypos, cos(2 * pi * (1:nelem(idx)) / nelem(idx)) * radius(idx)];
  ok = false;
  count = 0;
  while ~ok
    count = count + 1;
    phirand = (2 * pi / nelem(idx)) * (rand(1, nelem(idx)) - 0.5);
    rrand = 0.1 * radius(idx) * (rand(1, nelem(idx)) - 0.5);
    sel = sum([6. nelem(1:idx-1)]) + 1:sum([6, nelem(1:idx)]);
    r0 = sqrt(xpos(sel).^2 + ypos(sel).^2);
    phi0 = atan2(ypos(sel), xpos(sel));
    xposr = (r0 + rrand) .* sin(phi0 + phirand);
    yposr = (r0 + rrand) .* cos(phi0 + phirand);
    u = meshgrid([xpos(1:sel(1)-1), xposr]) - meshgrid([xpos(1:sel(1)-1), xposr]).';
    v = meshgrid([ypos(1:sel(1)-1), yposr]) - meshgrid([ypos(1:sel(1)-1), yposr]).';
    dist = sqrt(u(:).^2 + v(:).^2);
    if (sum(dist(dist > 0) < 2.4) == 0 | count >= 10 )
      disp(num2str(count));
      if (count < 10)
        xpos(sel) = xposr;
        ypos(sel) = yposr;
      end
      ok = true;
    end
  end
end

% scenario JDB
dmin = 2.4;
dmax = 20;
rmin = dmin;
rmax = 35;
nelem = 90;
incr = (rmax - rmin) / nelem;
dphimin = asin(dmin / rmin);
dphimax = asin(dmax / rmax);
xpos2 = (rmin + incr * (1:nelem-1)) .* sin(cumsum(dphimin:(dphimax - dphimin)/(nelem-2):dphimax));
ypos2 = (rmin + incr * (1:nelem-1)) .* cos(cumsum(dphimin:(dphimax - dphimin)/(nelem-2):dphimax));

