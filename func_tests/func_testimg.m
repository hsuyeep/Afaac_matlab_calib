% Script to compare the imaging between direct and FFT imaging.
% pep/06May15

close all; clear all;
run ~/startup.m
fid = fopen ('/dop312_0/prasad/dop147_0_AARTFAAC/20Nov13/r03/SB002_LBA_OUTER_8b2sbr03_1ch_1850-2127.bin', 'rb');
[acc, tobs, fobs] = readms2float (fid, -1, -1, 288);
load ('poslocal_outer.mat', 'posITRF', 'poslocal');

flagant = [129, 140, 149];
acc_wh = acc ./sqrt(diag(acc) * diag(acc)');
mask = ones (288);
mask (flagant,:) = 0;
mask (:,flagant) = 0;
acfl_wh = reshape(acc_wh(mask(:) == 1), [285, 285]);
l = [-1:0.01:1];

rem_ants = setdiff ([1:288], flagant);
posITRF_fl (:,1) = posITRF(rem_ants, 1);
posITRF_fl (:,2) = posITRF(rem_ants, 2);

% DFT map
map = acm2skyimage (conj(acc), posITRF(:,1), posITRF(:,2), fobs, l, l);
imagesc(l,l,fliplr(abs(map))); colorbar;
set(gca, 'FontSize', 16);
title (sprintf ('Uncalib. DFT. Time: %s, Freq: %.2f', datestr(mjdsec2datenum(tobs)), fobs));
axis equal
axis tight
set (gca, 'YDir', 'Normal'); % To match orientation with station images
set (gca, 'XDir', 'Reverse'); % To match orientation with station images
ylabel('South $\leftarrow$ m $\rightarrow$ North', 'interpreter', 'latex');
xlabel('East $\leftarrow$ l $\rightarrow$ West', 'interpreter', 'latex');
print ([mfilename '_dft.png'], '-dpng', '-r300');
fprintf (1, 'uncalib covariance DFT:  pixel range: %d-%d\n', max(max(abs(map))), min(min(abs(map))));

map_wh = acm2skyimage (conj(acc_wh), posITRF(:,1), posITRF(:,2), fobs, l, l);
imagesc(l,l,fliplr(abs(map_wh))); colorbar;
set(gca, 'FontSize', 16);
title (sprintf ('Uncalib. DFT whiten. Time: %s, Freq: %.2f', datestr(mjdsec2datenum(tobs)), fobs));
axis equal
axis tight
set (gca, 'YDir', 'Normal'); % To match orientation with station images
set (gca, 'XDir', 'Reverse'); % To match orientation with station images
ylabel('South $\leftarrow$ m $\rightarrow$ North', 'interpreter', 'latex');
xlabel('East $\leftarrow$ l $\rightarrow$ West', 'interpreter', 'latex');
print ([mfilename '_dft_whiten.png'], '-dpng', '-r300');
fprintf (1, 'uncalib corr.coef (whiten)DFT:  pixel range: %d-%d\n', max(max(abs(map_wh))), min(min(abs(map_wh))));

map_wh_fl = acm2skyimage (conj(acfl_wh), posITRF_fl(:,1), posITRF_fl(:,2), fobs, l, l);
imagesc(l,l,fliplr(abs(map_wh_fl))); colorbar;
set(gca, 'FontSize', 16);
title (sprintf ('Uncalib. DFT whiten flag. Time: %s, Freq: %.2f', datestr(mjdsec2datenum(tobs)), fobs));
axis equal
axis tight
set (gca, 'YDir', 'Normal'); % To match orientation with station images
set (gca, 'XDir', 'Reverse'); % To match orientation with station images
ylabel('South $\leftarrow$ m $\rightarrow$ North', 'interpreter', 'latex');
xlabel('East $\leftarrow$ l $\rightarrow$ West', 'interpreter', 'latex');
print ([mfilename '_dft_whiten_flag.png'], '-dpng', '-r300');
fprintf (1, 'uncalib corr.coef (whiten+flag)DFT:  pixel range: %d-%d\n', max(max(abs(map_wh_fl))), min(min(abs(map_wh_fl))));

% FFT map
gparm.type = 'pillbox';
gparm.lim  = 0;
gparm.duv = 0.5;                % Default, reassigned from fobs. of obs. to
                                % image just the full Fov (-1<l<1)
gparm.Nuv = 500;                % size of gridded visibility matrix
gparm.uvpad = 512;              % specifies if any padding needs to be added
gparm.fft  = 1;

uloc = meshgrid (poslocal(:,1)) - meshgrid (poslocal(:,1)).';
vloc = meshgrid (poslocal(:,2)) - meshgrid (poslocal(:,2)).';
[ra, fft_1, calvis, l, m] = fft_imager_sjw_radec (acc(:), uloc(:), vloc(:), gparm, [], [], tobs, fobs, 0);
imagesc(l,m,real(fft_1)); colorbar;
set(gca, 'FontSize', 16);
title (sprintf ('Uncalib. FFT Time: %s, Freq: %.2f', datestr(mjdsec2datenum(tobs)), fobs));
axis equal
axis tight
set (gca, 'YDir', 'Normal'); % To match orientation with station images
set (gca, 'XDir', 'Reverse'); % To match orientation with station images
ylabel('South $\leftarrow$ m $\rightarrow$ North', 'interpreter', 'latex');
xlabel('East $\leftarrow$ l $\rightarrow$ West', 'interpreter', 'latex');
print ([mfilename '_fft.png'], '-dpng', '-r300');
fprintf (1, 'uncalib covariance FFT:  pixel range: %d-%d\n', max(max(real(fft_1))), min(min(real(fft_1))));

[ra, fft_wh, calvis, l, m] = fft_imager_sjw_radec (acc_wh(:), uloc(:), vloc(:), gparm, [], [], tobs, fobs, 0);
imagesc(l,m,real(fft_wh)); colorbar;
set(gca, 'FontSize', 16);
title (sprintf ('Uncalib. FFT whiten Time: %s, Freq: %.2f', datestr(mjdsec2datenum(tobs)), fobs));
axis equal
axis tight
set (gca, 'YDir', 'Normal'); % To match orientation with station images
set (gca, 'XDir', 'Reverse'); % To match orientation with station images
ylabel('South $\leftarrow$ m $\rightarrow$ North', 'interpreter', 'latex');
xlabel('East $\leftarrow$ l $\rightarrow$ West', 'interpreter', 'latex');
print ([mfilename '_fft_whiten.png'], '-dpng', '-r300');
fprintf (1, 'uncalib corr.coef(whiten)FFT:  pixel range: %d-%d\n', max(max(real(fft_wh))), min(min(real(fft_wh))));


%%%%%%%%%%%%%%%% Calibrate %%%%%%%%%%%%%%%%%%%
sol = pelican_sunAteamsub (acc, tobs, fobs, eye(288), flagant, 0, 1, [], [], 'poslocal_outer.mat', [], []);
acccal = (sol.gainsol*sol.gainsol') .* acc;
[ra, fft_cal, calvis, l, m] = fft_imager_sjw_radec (acccal(:), uloc(:), vloc(:), gparm, [], [], tobs, fobs, 0);
imagesc(l,m,real(fft_cal)); colorbar;
set(gca, 'FontSize', 16);
title (sprintf ('Calib. FFT Time: %s, Freq: %.2f', datestr(mjdsec2datenum(tobs)), fobs));
axis equal
axis tight
set (gca, 'YDir', 'Normal'); % To match orientation with station images
set (gca, 'XDir', 'Reverse'); % To match orientation with station images
ylabel('South $\leftarrow$ m $\rightarrow$ North', 'interpreter', 'latex');
xlabel('East $\leftarrow$ l $\rightarrow$ West', 'interpreter', 'latex');
print ([mfilename '_cal_fft.png'], '-dpng', '-r300');
fprintf (1, 'Calib covariance FFT:  pixel range: %d-%d\n', max(max(real(fft_cal))), min(min(real(fft_cal))));

[uloc_fl, vloc_fl] = gen_flagged_uvloc (uloc, vloc, flagant); 
acccalsign = (sol.gainsol(rem_ants)*sol.gainsol(rem_ants)') .* acfl_wh - sol.sigman;
[ra, fft_calsign, calvis, l, m] = fft_imager_sjw_radec (acccalsign(:), uloc_fl(:), vloc_fl(:), gparm, [], [], tobs, fobs, 0);
imagesc(l,m,real(fft_calsign)); colorbar;
set(gca, 'FontSize', 16);
title (sprintf ('Calib. FFT, noise sub Time: %s, Freq: %.2f', datestr(mjdsec2datenum(tobs)), fobs));
axis equal
axis tight
set (gca, 'YDir', 'Normal'); % To match orientation with station images
set (gca, 'XDir', 'Reverse'); % To match orientation with station images
ylabel('South $\leftarrow$ m $\rightarrow$ North', 'interpreter', 'latex');
xlabel('East $\leftarrow$ l $\rightarrow$ West', 'interpreter', 'latex');
print ([mfilename '_cal_fft_nonoise.png'], '-dpng', '-r300');
fprintf (1, 'Calib covariance, noise sub FFT:  pixel range: %d-%d\n', max(max(real(fft_calsign))), min(min(real(fft_calsign))));

[ra, fft_calnoa, calvis, l, m] = fft_imager_sjw_radec (sol.calvis+diag(diag(sol.sigman)), uloc_fl(:), vloc_fl(:), gparm, [], [], tobs, fobs, 0);
imagesc(l,m,real(fft_calnoa)); colorbar;
set(gca, 'FontSize', 16);
title (sprintf ('Calib. FFT, noAteam Time: %s, Freq: %.2f', datestr(mjdsec2datenum(tobs)), fobs));
axis equal
axis tight
set (gca, 'YDir', 'Normal'); % To match orientation with station images
set (gca, 'XDir', 'Reverse'); % To match orientation with station images
ylabel('South $\leftarrow$ m $\rightarrow$ North', 'interpreter', 'latex');
xlabel('East $\leftarrow$ l $\rightarrow$ West', 'interpreter', 'latex');
print ([mfilename '_cal_fft_noa.png'], '-dpng', '-r300');
fprintf (1, 'Calib covariance, noise+Ateamsub  FFT:  pixel range: %d-%d\n', max(max(real(fft_calnoa))), min(min(real(fft_calnoa))));
