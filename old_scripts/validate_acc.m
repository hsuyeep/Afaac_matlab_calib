function validate_acc(acc, calfile, posfile, rcumode, sbsel, tobs)

% validate_acc(acc, calfile, posfile, rcumode, sbsel, tobs)
%
% Validation of the calibration table using an array covariance cube
% provided by the user (either obtained from a calibration observation or
% composed from multiple xst-files) by comparison of images before and
% after calibration. The images produced are averaged over the use
% specified subband selection.
%
% Arguemens
% acc     : Nrcu x Nrcu x Nsb array covariance cube
% calfile : calibration table file as produced by writeCalTable
% posfile : AntennaField.conf file of appropriate station
% rcumode : RCU mode number
% sbsel   : Nsb x 1 vector with subband selection
% t_obs   : time of observation as dete number
%
% SJW, October 2011

% read the data
Nsb = size(acc, 3);
Nrcu = size(acc, 1);
if (Nrcu == 192)
    EU = true;
else
    EU = false;
end

% read the calibration table
[calx, caly] = readCalTable(calfile);
cal = zeros(Nrcu, 512);
cal(1:2:end, :) = calx;
cal(2:2:end, :) = caly;

% determine subband center frequencies
if (rcumode == 1 || rcumode == 2 || rcumode == 3 || rcumode == 4)
    freq = 0:1e8/512:1e8-1e-3;
elseif (rcumode == 5)
    freq = 1e8:1e8/512:2e8-1e-3;
elseif (rcumode == 6)
    freq = 160e6:80e6/512:240e6-1e-3;
elseif (rcumode == 7)
    freq = 2e8:1e8/512:3e8-1e-3;
else
    disp(['RCU mode numer ' num2str(rcumode) ' is not a valid RCU mode.']);
    return
end
freq = freq(sbsel);

% calibrate acc
acccal = zeros(size(acc));
for sb = 1:Nsb
    acccal(:, :, sb) = conj(cal(:, sbsel(sb)) * cal(:, sbsel(sb))') .* acc(:, :, sb);
end

% read antenna positions
[pos, ~, ~, ~] = parseAntennaField(posfile, rcumode, EU);
xpos = pos(:, 1);
ypos = pos(:, 2);

% produce sky maps
l = -1:0.02:1;
m = -1:0.02:1;
skymapx = acm2skyimage(acc(1:2:end, 1:2:end, :), xpos, ypos, freq, l, m);
skymapxcal = acm2skyimage(acccal(1:2:end, 1:2:end, :), xpos, ypos, freq, l, m);
skymapy = acm2skyimage(acc(2:2:end, 2:2:end, :), xpos, ypos, freq, l, m);
skymapycal = acm2skyimage(acccal(2:2:end, 2:2:end, :), xpos, ypos, freq, l, m);

% show the result
dist = sqrt(meshgrid(l).^2 + meshgrid(m).'.^2);
mask = NaN(size(dist));
mask(dist < 1) = 1;
figure
set(gcf, 'Position', [0, 500, 1120 840]);
subplot(2, 2, 1);
imagesc(l, m, mean(skymapx, 3).' .* mask);
set(gca, 'FontSize', 16, 'YDir', 'normal', 'XDir', 'reverse');
set(colorbar, 'FontSize', 16);
xlabel('East \leftarrow l \rightarrow West');
ylabel('South \leftarrow m \rightarrow North');
title('Uncalibrated map in x-polarization');
subplot(2, 2, 2);
imagesc(l, m, mean(skymapy, 3).' .* mask);
set(gca, 'FontSize', 16, 'YDir', 'normal', 'XDir', 'reverse');
set(colorbar, 'FontSize', 16);
xlabel('East \leftarrow l \rightarrow West');
ylabel('South \leftarrow m \rightarrow North');
title('Uncalibrated map in y-polarization');
subplot(2, 2, 3);
imagesc(l, m, mean(skymapxcal, 3).' .* mask);
set(gca, 'FontSize', 16, 'YDir', 'normal', 'XDir', 'reverse');
set(colorbar, 'FontSize', 16);
xlabel('East \leftarrow l \rightarrow West');
ylabel('South \leftarrow m \rightarrow North');
title('Calibrated map in x-polarization');
subplot(2, 2, 4);
imagesc(l, m, mean(skymapycal, 3).' .* mask);
set(gca, 'FontSize', 16, 'YDir', 'normal', 'XDir', 'reverse');
set(colorbar, 'FontSize', 16);
xlabel('East \leftarrow l \rightarrow West');
ylabel('South \leftarrow m \rightarrow North');
title('Calibrated map in y-polarization');
h = text(1.4, 4.15, ['Observed in RCU mode ' num2str(rcumode) ' at ' datestr(tobs) ' UTC']);
set(h, 'FontSize', 16, 'HorizontalAlignment', 'center');
