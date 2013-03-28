function [cal, freq, skymap] = CS10pipeline(acc, t_obs, flags, N, xpos, ypos, l, m, freq, elembeam, lbeam, mbeam)

% RFI detection
chsel = RFIdetection(acc, N);
freq = freq(chsel);
acc = acc(:, :, chsel);

% flagging
for idx = 1:length(freq);
    acc(:, :, idx) = acc(:, :, idx) .* (1 - flags);
end

% calibration
c = 2.99792e8;
restriction = 4; % baseline restriction in wavelengths
lat = 52 + 54/60 + 42.6/3600; % geographical location CS010
lon = 6 + 52/60 + 3.17/3600;
cal = zeros(length(freq), length(xpos));
for idx = 1:length(freq)
    disp([num2str(idx) ' of ' num2str(length(freq))]);
    [lsrc, msrc, flux] = makesrcmodel(t_obs(idx), lon, lat);
    att = interp2(lbeam, mbeam, squeeze(elembeam(:, :, idx)), lsrc, msrc);
    appflux = flux .* att;
    srcmat = [appflux, lsrc, msrc];
    % note: baseline restriction limited to 24 meters at lowest wavelengths
    % to avoid NaN in computation of alpha.
    [cal(idx, :), quality] = cal_msrc(acc(:, :, idx), freq(idx), xpos, ypos, min([restriction, 20 / (c / freq(idx))]), srcmat);
end

% apply calibration
acmcal = zeros(size(acc));
for idx = 1:length(freq)
    acmcal(:, :, idx) = (cal(idx, :)' * cal(idx, :)) .* squeeze(acc(:, :, idx));
end

%for idx = 1:length(freq)
%    [lsrc, msrc, flux] = makesrcmodel(t_obs(idx), lon, lat);
%    att = interp2(lbeam, mbeam, squeeze(elembeam(:, :, idx)), lsrc, msrc);
%    appflux = flux .* att;
%    srcmat = [appflux, lsrc, msrc];
%    [caltest(idx, :), quality] = cal_msrc(acmcal(:, :, idx), freq(idx), xpos, ypos, min(restriction, 24 / (c / freq(idx))), srcmat);
%end

% imaging
skymap = acm2skyimage(acmcal, xpos, ypos, freq, l, m);
