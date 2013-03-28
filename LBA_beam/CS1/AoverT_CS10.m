% Aeff/Tsys computations for array of x-dipoles on CS10
% using failure mode research done on March 14

%snapshot = 79;
c = 2.99792e8;
restriction = min(25 * freq / c, 4);
l = -1:0.005:1;
m = -1:0.005:1;

% element beam pattern
load cs10_fieldi.mat
[lgrid, mgrid] = meshgrid(l, m);
dist = sqrt(lgrid.^2 + mgrid.^2);
thetai = asin(dist);
thetai(dist >= 1) = 0;
phii = mod(atan2(lgrid, mgrid), 2 * pi);
E_11 = zeros(length(sel), size(FieldI(1).PHI, 1), size(FieldI(1).PHI, 2), length(FieldI(1).Freq));
for fidx = 1:length(FieldI);
  for idx = 1:length(sel);
    J_11 = FieldI(fidx).E(:, :, 2, selx(sel(idx)));
    J_12 = FieldI(fidx).E(:, :, 3, selx(sel(idx)));
    E_11(idx, :, :, fidx) = 0.5 * (J_11 .* conj(J_11) + J_12 .* conj(J_12));
  end
    simfreq(fidx) = FieldI(fidx).Freq;
end
Eavg11 = squeeze(mean(E_11, 1));
elembeam_xx = zeros(length(l), length(m), length(freq));
for idx = 1:length(freq)
  elembeam_xx(:, :, idx) = interp3(FieldI(1).PHI(1, :), FieldI(1).THETA(:, 1), simfreq, Eavg11, phii, thetai, freq(idx) * ones(size(phii)));
end

% calibration
load srclist3C.mat
[lsrc, msrc] = radectolm(srclist3C.alpha, srclist3C.delta, JulianDay(t_obs(snapshot)));
up = ~isnan(lsrc);
lsrc = lsrc(up);
msrc = msrc(up);
cal0 = zeros(length(freq), length(sel));
quality = zeros(length(freq));
for idx = 1:length(freq)
  att = interp2(l, m, squeeze(elembeam_xx(:, :, idx)), lsrc, msrc);
  Isrc_app = srclist3C.flux(up) .* att;
  srcmat = [Isrc_app, lsrc, msrc];
  [cal0(idx, :), quality(idx)] = cal_msrc(data{snapshot}(selx(sel), selx(sel), idx), freq(idx), xpos(CS10idx(sel)), ypos(CS10idx(sel)), restriction(idx), srcmat);
end
cal0 = cal0 ./ abs(cal0);

acmcal = zeros(size(data{snapshot}(selx(sel), selx(sel), :)));
for idx = 1:length(freq)
  acmcal(:, :, idx) = (cal0(idx, :)' * cal0(idx, :)) .* squeeze(data{snapshot}(selx(sel), selx(sel), idx));
end

%verification
%[cal1, quality] = calacm(acmcal, freq, xpos(CS10idx(sel)), ypos(CS10idx(sel)), restriction, t_obs(snapshot));

% imaging
[lCas, mCas] = radectolm(srclist3C.alpha(324), srclist3C.delta(324), JulianDay(t_obs(snapshot)));
baseline = sqrt((meshgrid(xpos(CS10idx(sel))) - meshgrid(xpos(CS10idx(sel))).').^2 + (meshgrid(ypos(CS10idx(sel))) - meshgrid(ypos(CS10idx(sel))).').^2);
limg = l(find(abs(l - lCas) < 0.1));
mimg = m(find(abs(m - mCas) < 0.1));
skymapxcal = zeros(length(limg), length(mimg), length(freq));
for idx = 1:length(freq)
  skymapxcal(:, :, idx) = acm2skyimage(squeeze(acmcal(:, :, idx)) .* (baseline > restriction(idx) * c / freq(idx)), xpos(CS10idx(sel)), ypos(CS10idx(sel)), freq(idx), limg, mimg);
end

for idx = 1:length(freq)
  PCas(idx) = max(max(skymapxcal(:, :, idx))) / sum(sum(baseline > restriction(idx) * c / freq(idx)));
  acpower = diag(squeeze(acmcal(:, :, idx)));
  Pac(idx) = mean(acpower);
end

% flux Cas A from Baars et al.
logS1965 = 5.625 - 0.634 * log10(freq / 1e6) - 0.023 * log10(freq / 1e6).^2;
S1965 = 10.^logS1965;
dflux = 0.97 - 0.30 * log10(freq / 1e9);
SCas = S1965 .* (1 - dflux / 100).^(2007-1965) * 1e-26;
mask = ones(length(l), length(m));
mask(dist >= 1) = 0;
for idx = 1:length(freq)
  gCas(idx) = interp2(lgrid, mgrid, squeeze(elembeam_xx(:, :, idx)), lCas, mCas);
  gmax(idx) = max(max(squeeze(elembeam_xx(:, :, idx)) .* mask));
end
SCas = (gCas ./ gmax) .* SCas / sqrt(1 - lCas^2 - mCas^2);

% compute Ae as funtion of frequency
weight = sqrt(1 - meshgrid(l).^2 - meshgrid(m).'.^2);
weight(weight == 0) = 1e-10;
sum(sum(mask ./ weight)) / sum(sum(mask))
Omega_e = zeros(size(freq));
for idx = 1:length(freq)
  Omega_e(idx) = sum(sum(squeeze(elembeam_xx(:, :, idx) / gmax(idx)) .* mask ./ weight)) / sum(sum(mask)) * pi;
end
lambda = c ./ freq;
Ae = lambda.^2 ./ Omega_e;

% Tsys
k = 1.38e-23;
Tsys = Ae .* (SCas / (2 * k)) .* ((Pac - PCas) ./ PCas);
