freqidx = 370:380;
for idx = 1:length(freqidx);
   acm = acc(:, :, freqidx(idx));
   acmx = acm;
   acmxac = diag(acmx);
   acmxcc = triu(acmx, 1);
   acmx = acmxcc' + diag(acmxac) + acmxcc;
   acmxfull(:, :, idx) = acmx(2:2:end, 2:2:end);
end

sel = [0:19, 22:33, 36:37, 40:43, 46:61, 64:93] + 1;
sel = sel(2:2:end)/2;

l = -1:0.02:1;
m = -1:0.02:1;

skymap0 = acm2skyimage(acmxfull(sel, sel, :), xpos(CS10idx(sel)), ypos(CS10idx(sel)), freq(freqidx), l, m);
dist = sqrt(meshgrid(l, m).^2 + meshgrid(m, l).'.^2);
mask = dist < 1;
mask = mask .* 1.0;
mask(mask == 0) = NaN;
pcolor(m, l, squeeze(sum(skymap0, 3)) .* mask.');
shading flat