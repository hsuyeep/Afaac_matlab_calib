% Script to compare LBA_INNER and LBA_OUTER PSFs with and without optimization.

flagant=[47:48:288, 48:48:288];
freq=60000000;
basedir = '~/Documents/AARTFAAC_Project_SW_system_plan/afaac_psfshape/';
tparm=[]; wparm=[]; gparm=[]; deb=0;
[l,m,psf_in,w_in, intap_in, outtap_in] = genarraypsf ('poslocal_inner.mat', flagant, freq, tparm, wparm, gparm, deb);
psf_in = psf_in/(max(max(psf_in)));
mask = NaN (size (psf_in));
mask (meshgrid(l).^2 + meshgrid(m).'.^2 < 1) = 1;

% Plot antenna positions
load ('poslocal_inner.mat', 'poslocal');
central_ant = [poslocal(1,1), poslocal(1,2)];
plot (poslocal(:,1)-central_ant(1), poslocal(:,2)-central_ant(2), '.');
hold on;
clear poslocal;
load ('poslocal_outer.mat', 'poslocal');
plot (poslocal(:,1)-central_ant(1), poslocal(:,2)-central_ant(2), 'r.');
legend ('LBA\_INNER', 'LBA\_OUTER');
xlabel ('X (m)'); ylabel ('Y (m)');
grid on; axis tight;
print ([basedir, 'LBA_INNER_outer_antpos.eps'], '-depsc');

clf;
imagesc (l,m,10*log10(psf_in.*mask));
caxis ([-50, 0]);
colorbar;
set(gca, 'fontsize', 16);
title ('LBA\_INNER PSF, 60MHz');
axis equal
axis tight
set (gca, 'ydir', 'normal'); % to match orientation with station images
set (gca, 'xdir', 'reverse'); % to match orientation with station images
% ylabel('South $\leftarrow$ m $\rightarrow$ North', 'interpreter', 'latex');
% xlabel('East $\leftarrow$ l $\rightarrow$ West', 'interpreter', 'latex');
ylabel('South \leftarrow m \rightarrow North');
xlabel('East \leftarrow l \rightarrow West');
print ([basedir, 'LBA_INNER_PSF.eps'], '-depsc');

% PSF profile
clf;
subplot (2,1,1);
plot (l, 10*log10(psf_in (:,128)));
ylim ([-50 0]);
xlabel ('l'); ylabel ('dB');
grid on; axis tight;
text (0.3, -5, '\bf LBA\_INNER, l');
% print ([basedir, 'LBA_INNER_PSF_lprofile.eps'], '-depsc');

subplot (2,1,2);
plot (m, 10*log10(psf_in (128,:)));
ylim ([-50 0]);
grid on; axis tight;
text (0.3, -5, '\bf LBA\_INNER, m');
xlabel ('l/m'); ylabel ('dB');
samexaxis ('join');
print ([basedir, 'LBA_INNER_PSF_profile.eps'], '-depsc');


flagant=[];
[l,m,psf_out,w_out, outtap_out, outtap_out] = genarraypsf ('poslocal_outer.mat', flagant, freq, tparm, wparm, gparm, deb);
psf_out = psf_out/(max(max(psf_out)));
figure;

% Plot antenna positions
load ('poslocal_outer.mat', 'poslocal');
plot (poslocal(:,1)-poslocal(1,1), poslocal(:,2)-poslocal(1,2), '.');
grid on; axis tight;
xlabel ('X (m)'); ylabel ('Y (m)');
print ([basedir, 'LBA_OUTER_antpos.eps'], '-depsc');

clf;
imagesc (l,m,10*log10(psf_out.*mask));
caxis ([-50 0]);
set(gca, 'fontsize', 16);
title ('LBA\_OUTER PSF, 60MHz');
axis equal
axis tight
set (gca, 'ydir', 'normal'); % to match orientation with station images
set (gca, 'xdir', 'reverse'); % to match orientation with station images
% ylabel('south $\leftarrow$ m $\rightarrow$ north', 'interpreter', 'latex');
% xlabel('east $\leftarrow$ l $\rightarrow$ west', 'interpreter', 'latex');
ylabel('South \leftarrow m \rightarrow North');
xlabel('East \leftarrow l \rightarrow West');
colorbar;
print ([basedir, 'LBA_OUTER_PSF.eps'], '-depsc');

% PSF profile
clf;
subplot (211);
plot (l, 10*log10(psf_out (:,128)));
ylim ([-50 0]);
xlabel ('l'); ylabel ('dB');
grid on; axis tight;
text (0.3, -5, '\bf LBA\_OUTER, l');
% print ([basedir, 'LBA_OUTER_PSF_lprofile.eps'], '-depsc');

subplot (212);
plot (m, 10*log10(psf_out (128,:)));
ylim ([-50 0]);
grid on; axis tight;
text (0.3, -5, '\bf LBA\_OUTER, m');
xlabel ('m'); ylabel ('dB');
samexaxis ('join');
print ([basedir, 'LBA_OUTER_PSF_profile.eps'], '-depsc');

%%%%%%%%%%%%%%%%%% Optimized PSF %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generate the optimized PSF after applying weighting, taper, GCF.
% 1. LBA_INNER uniform weighting (10lambda), outer taper
% Generate taper parameters
tparm.type = 'Gaussian';
tparm.minlambda = 10;    % Spatial freq. cutoff lower limit.NOTE: Units of lambda.
tparm.maxmeters = 350;	 % Spatial freq. cutoff upper limit. NOTE:Units of meters.
tparm.pa(1) = 10;	 % Inner taper gaussian sigx. NOTE: Units of lambda. 
tparm.pa(2) = tparm.pa(1);% Inner taper gaussian sigy
tparm.pa(3) = 20; 	     % Outer taper Gaussian sigx. NOTE: Units of lambda.
tparm.pa(4) = tparm.pa(3); % Outer taper Gaussian sigy.

% Generate weight parameters
wparm.type = 'uniform'; wparm.cellrad = 10; % lambda

% Generate parameters for a gaussian GCF
gaugcf.type = 'Gaussian';
gaugcf.duv = 0.5; % Wavelength
gaugcf.Nuv = 240; % Points
gaugcf.uvpad = 256; 
gaugcf.lim = 1.5; % Limit of gaussian kernel in wavelength units.
gaugcf.pa(1) = 0.5; % 2-D gaussian x/ysigma, in wavelength units.
gaugcf.pa(2) = 0.5;
gaugcf.fft = 1;

flagant=[47:48:288, 48:48:288];
[l,m,psf_in_opt,weight, intap, outtap, uvdist] = genarraypsf ('poslocal_inner.mat', flagant, freq, tparm, wparm, gaugcf, deb);

clf;
psf_in_opt = psf_in_opt/(max(max(psf_in_opt)));
imagesc (l,m,10*log10(abs(psf_in_opt.*mask)));
caxis ([-50 0]);
colorbar;
set(gca, 'fontsize', 16);
title ('LBA\_INNER PSF, 60MHz');
axis equal
axis tight
set (gca, 'ydir', 'normal'); % to match orientation with station images
set (gca, 'xdir', 'reverse'); % to match orientation with station images
% ylabel('South $\leftarrow$ m $\rightarrow$ North', 'interpreter', 'latex');
% xlabel('East $\leftarrow$ l $\rightarrow$ West', 'interpreter', 'latex');
ylabel('South \leftarrow m \rightarrow North');
xlabel('East \leftarrow l \rightarrow West');
print ([basedir, 'LBA_INNER_PSF_opt.eps'], '-depsc');

% PSF profile
clf;
subplot (2,1,1);
plot (l, 10*log10(psf_in_opt (:,128)));
ylim ([-50 0]);
xlabel ('l'); ylabel ('dB');
grid on; axis tight;
text (0.3, -5, '\bf LBA\_INNER, l');
% print ([basedir, 'LBA_INNER_PSF_lprofile.eps'], '-depsc');

subplot (2,1,2);
plot (m, 10*log10(psf_in_opt (128,:)));
ylim ([-50 0]);
grid on; axis tight;
text (0.3, -5, '\bf LBA\_INNER, m');
xlabel ('l/m'); ylabel ('dB');
samexaxis ('join');
print ([basedir, 'LBA_INNER_PSF_opt_profile.eps'], '-depsc');

% Optimization weights
subplot (211);
plot (uvdist, intap .* outtap, '.');
grid on; axis tight;
xlabel ('UV dist (m)');
ylabel ('Taper');
subplot (212);
plot (uvdist, weight, '.');
grid on; axis tight;
xlabel ('UV dist (m)');
ylabel ('Baseline density');
samexaxis ('join');
print ([basedir, 'LBA_INNER_PSF_opt_taper_weight.eps'], '-depsc');

%%%%%%%%%%%%%%%  Optimized PSF. LBA_OUTER %%%%%%%%%%%%%%%%%%%%%%
flagant=[];
[l,m,psf_out_opt, weight_out, intap_out, outtap_out, uvdist_out] = genarraypsf ('poslocal_outer.mat', flagant, freq, tparm, wparm, [], deb);

clf;
psf_out_opt = psf_out_opt/(max(max(psf_out_opt)));
imagesc (l,m,10*log10(abs(psf_out_opt.*mask)));
caxis ([-50 0]);
colorbar;
set(gca, 'fontsize', 16);
title ('LBA\_OUTER PSF, 60MHz');
axis equal
axis tight
set (gca, 'ydir', 'normal'); % to match orientation with station images
set (gca, 'xdir', 'reverse'); % to match orientation with station images
% ylabel('South $\leftarrow$ m $\rightarrow$ North', 'interpreter', 'latex');
% xlabel('East $\leftarrow$ l $\rightarrow$ West', 'interpreter', 'latex');
ylabel('South \leftarrow m \rightarrow North');
xlabel('East \leftarrow l \rightarrow West');
print ([basedir, 'LBA_OUTER_PSF_opt.eps'], '-depsc');

% PSF profile
clf;
subplot (2,1,1);
plot (l, 10*log10(psf_out_opt (:,128)));
ylim ([-50 0]);
xlabel ('l'); ylabel ('dB');
grid on; axis tight;
text (0.3, -5, '\bf LBA\_OUTER, l');
% print ([basedir, 'LBA_OUTER_PSF_lprofile.eps'], '-depsc');

subplot (2,1,2);
plot (m, 10*log10(psf_out_opt (128,:)));
ylim ([-50 0]);
grid on; axis tight;
text (0.3, -5, '\bf LBA\_OUTER, m');
xlabel ('l/m'); ylabel ('dB');
samexaxis ('join');
print ([basedir, 'LBA_OUTER_PSF_opt_profile.eps'], '-depsc');

%%%%%%%%%%%%%%%%%%%%%%%% Real data imaging based on optimized PSF %%%%%%%%%%%%%%%%%%%
% Data locations:
fid_inner = fopen ('~/WORK/AARTFAAC/Reobs/11Jul12/LBA_INNER_BAND60/SB000_LBA_INNER_BAND60_1ch.bin', 'rb');  
% Flag for this observation.
flagant = [1:12, 47:48:288, 48:48:288, 49:48:288]; 
[acc, tobs, freq] = readms2float (fid_inner, -1, -1, 288);
load ('poslocal_inner.mat', 'poslocal', 'posITRF');
antmask = zeros (size (acc));
posmask = zeros (size (poslocal));
rem_ants = length(acc) - length(flagant);
for ind = 1:length(flagant)
  antmask (flagant(ind), :) = 1; antmask (:,flagant(ind)) = 1;
  posmask (flagant(ind), :) = 1;
end
acc_fl = reshape (acc(antmask ~= 1), [rem_ants, rem_ants]);
% poslocal_fl = reshape(poslocal(posmask ~=1), [rem_ants, 3]); 
% posITRF_fl = reshape(posITRF(posmask ~=1), [rem_ants, 3]); 
[l,m,psf_in_opt,weight, intap, outtap, uvdist] = genarraypsf ('poslocal_inner.mat', flagant, freq, tparm, wparm, [], deb);

% apply optmization before carrying out calibration.
% NOTE: Since the calibration routine accepts an ACM with 288 ants, we create 
% such an ACM to satisfy it. Also, having the intap results in calibration failure
% due to a badly conditioned matrix, so we don't do that pre-calibration.
acc_opt = acc_fl (:) .* outtap ./weight;
acm_dummy = zeros (288);
acm_dummy(antmask == 0) = acc_opt;
sol = pelican_sunAteamsub (acm_dummy, tobs, fobs, eye(size (acm_dummy)), ... 
				flagant, 0, 1, [], [], 'poslocal_inner.mat', [], []);
l = [-1:0.01:1];
innermap = acm2skyimage (sol.calvis, poslocal_fl(:,1), poslocal_fl(:,2), fobs, l, l);


fid_outer = fopen ('~/WORK/AARTFAAC/Reobs/11Jul12/LBA_INNER_BAND60/SB000_LBA_INNER_BAND60_1ch.bin', 'rb');  
