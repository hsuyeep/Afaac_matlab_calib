% Script to attempt imaging from LBA_INNER/LBA_SPARSE configuration based on
% data acquired earlier.
% pep/28Apr15

close all; clear all;
array = 'inner';

basedir = '~/WORK/AARTFAAC/Afaac_matlab_calib/func_tests/';
% For imaging
load ('poslocal_inner.mat', 'poslocal', 'posITRF');

% Data locations:
fid_inner = fopen ('~/WORK/AARTFAAC/Reobs/11Jul12/LBA_INNER_BAND60/SB000_LBA_INNER_BAND60_1ch.bin', 'rb');  
% Single station processing
% flagant = [1:48, 95:288];

% Flag for this observation.
% flagant = [1:12, 47:48:288, 48:48:288, 49:48:288]; 

% Flag only outrigger dipoles
flagant = [47:48:288, 48:48:288]; %, 49:48:288]; 
doprint = 1;

[acc, tobs, fobs] = readms2float (fid_inner, -1, -1, 288);

antmask = zeros (size (acc));
posmask = zeros (size (poslocal));
rem_ants = length(acc) - length(flagant);
for ind = 1:length(flagant)
  antmask (flagant(ind), :) = 1; antmask (:,flagant(ind)) = 1;
  posmask (flagant(ind), :) = 1;
end
acc_fl = reshape (acc(antmask ~= 1), [rem_ants, rem_ants]);
poslocal_fl = reshape(poslocal(posmask ~=1), [rem_ants, 3]); 
posITRF_fl = reshape(posITRF(posmask ~=1), [rem_ants, 3]); 

% Generate UV coverage for selected antennas.
uloc_fl = meshgrid (poslocal_fl(:,1)) - meshgrid (poslocal_fl(:,1)).';
vloc_fl = meshgrid (poslocal_fl(:,2)) - meshgrid (poslocal_fl(:,2)).';

plot (uloc_fl(:), vloc_fl(:), '.');
h = gca;
set (h, 'FontWeight', 'bold');
xlabel ('\bf U(m)'); ylabel ('\bf V(m)');
title ('\bf UV coverage for LBA\_INNER array configuration');
print ([basedir 'LBA_INNER_uvcov.eps'], '-depsc');

%%%%%%%%%%% Setup weighting, taper and GCF
% Generate visibility taper with different parameters.
tparm.type = 'Gaussian';
tparm.minlambda = 0;
tparm.maxmeters = 350;
tparm.pa = [-1, -1, 25, 25];
[intap, outtap, den, mask, uvdist] = taper (posITRF_fl, tparm, -1, fobs, 4);

% Generate bline histogram distribution
hist (uvdist, 50);
xlabel ('\bf UVdist'); ylabel ('\bf Baselines');
xlim ([0 300]); ylim ([0 5000]);
h = gca;
set (h, 'FontWeight', 'bold');
print ([basedir 'LBA_INNER_bline_distr.eps'], '-depsc');

% Generate visibility weights with different parameters
wparm.type = 'Uniform';
wparm.cellrad = 10;
weight = genvisweight (posITRF_fl, wparm, 0);

% Generate parameters for a gaussian GCF
gaugcf.type = 'Gaussian';
gaugcf.duv = 0.5;
gaugcf.Nuv = 240;
gaugcf.uvpad = 256;
gaugcf.lim = 1.5; % Limit of gaussian kernel in wavelength units.
gaugcf.pa(1) = 0.5; % 2-D gaussian x/ysigma, in wavelength units.
gaugcf.pa(2) = 0.5;
gaugcf.fft = 0;

uloc_fl = meshgrid (poslocal_fl(:,1)) - meshgrid (poslocal_fl(:,1)).';
vloc_fl = meshgrid (poslocal_fl(:,2)) - meshgrid (poslocal_fl(:,2)).';
uvdist_fl = sqrt (uloc_fl(:).^2 + vloc_fl(:).^2);

[gauvis, gauind, gaurng] = genvisgrid (acc_fl(:), uloc_fl(:), vloc_fl(:), gaugcf, fobs, 1);
[ugrid, vgrid] = meshgrid (gaurng);
uvdist_gr = sqrt (ugrid(:).^2 + vgrid(:).^2);

% Generate parameters for a pillbox GCF.
pillgcf.type = 'pillbox';
pillgcf.duv = 0.5;
pillgcf.Nuv = 240;
pillgcf.uvpad = 256;
pillgcf.lim = 0; % Limit of gaussian kernel in wavelength units.
pillgcf.fft = 1;
[pillvis, pillind, pillrng] = genvisgrid (acc_fl(:), uloc_fl(:), vloc_fl(:), pillgcf, fobs, 1);


% Generate weights directly from the uvdist histogram
edges = linspace (min(uvdist), max(uvdist), 100);
[val, bin] = histc(uvdist, edges);
histweights = 1./(val(bin).*bin);
count = accumarray(bin, histweights);
bar(count);

% Create uncalibrated image map
l = [-1:0.01:1];
innermap = acm2skyimage (acc_fl, poslocal_fl(:,1), poslocal_fl(:,2), fobs, l, l);
mask = NaN(size (innermap));
mask(meshgrid(l).^2 + meshgrid(l).'.^2 < 1) = 1;
imagesc (l, l, innermap.* mask);
colorbar;
set(gca, 'FontSize', 16);
title (sprintf ('LBA_INNER, uncalib Time: %s Freq: %.2f', datestr(mjdsec2datenum(tobs)), fobs));
axis equal
axis tight
set (gca, 'YDir', 'Normal'); % To match orientation with station images
set (gca, 'XDir', 'Reverse'); % To match orientation with station images

ylabel('South $\leftarrow$ m $\rightarrow$ North', 'interpreter', 'latex');
xlabel('East $\leftarrow$ l $\rightarrow$ West', 'interpreter', 'latex');
drawnow ();
if (doprint == 1)
    print ('~/Documents/AARTFAAC_Project_SW_system_plan/afaac_psfshape/LBA_INNER_analysis/afaac_uncal.png', '-dpng', '-r300');
end;
% Create uncalibrated single station maps
for stat = 1:6
    stflag = setdiff ([1:288], setdiff ((stat-1)*48 + [1:48], flagant));
    antmask = zeros (size (acc));
    posmask = zeros (size (poslocal));
    rem_ants = length(acc) - length(stflag);
    for ind = 1:length(stflag)
        antmask (stflag(ind), :) = 1; antmask (:,stflag(ind)) = 1;
        posmask (stflag(ind), :) = 1;
    end
    acc_st = reshape (acc(antmask ~= 1), [rem_ants, rem_ants]);
    poslocal_fl = reshape(poslocal(posmask ~=1), [rem_ants, 3]); 
    
    statmap = acm2skyimage (acc_st, poslocal_fl(:,1), poslocal_fl(:,2), fobs, l, l);
    imagesc (l, l, statmap.* mask);
    colorbar;
    set(gca, 'FontSize', 16);
    title (sprintf ('St. %d, uncalib Time: %s Freq: %.2f', stat, datestr(mjdsec2datenum(tobs)), fobs));
    axis equal
    axis tight
    set (gca, 'YDir', 'Normal'); % To match orientation with station images
    set (gca, 'XDir', 'Reverse'); % To match orientation with station images

    ylabel('South $\leftarrow$ m $\rightarrow$ North', 'interpreter', 'latex');
    xlabel('East $\leftarrow$ l $\rightarrow$ West', 'interpreter', 'latex');
    drawnow ();
    if (doprint == 1)
        print (sprintf ('~/Documents/AARTFAAC_Project_SW_system_plan/afaac_psfshape/LBA_INNER_analysis/st%d_uncal.png', stat), '-dpng', '-r300');
    end;
end;

% Calibrate ACM
sol = pelican_sunAteamsub (acc, tobs, fobs, eye(size (acc)), ... 
				flagant, 0, 1, [], [], 'poslocal_inner.mat', [], []);
clf
subplot (2,1,1);
plot (abs(sol.gainsol));
xlabel ('Ant'); ylabel ('Gain abs.');
subplot (2,1,2);
plot (angle(sol.gainsol));
xlabel ('Ant'); ylabel ('Gain phase (rad).');
drawnow ();
if (doprint == 1)
        print ('~/Documents/AARTFAAC_Project_SW_system_plan/afaac_psfshape/LBA_INNER_analysis/calsol.png', '-dpng', '-r300');
end;            
clf;

% First put the calibrated ACM into a 288 size ACM
acm_cal = zeros (288);
flag_mask = zeros (288);
flag_mask(flagant,:) = 1;
flag_mask(:,flagant) = 1;
 
% Adding noise estimation, else the single station calibrated images are zero.
% acm_cal (flag_mask == 0) = sol.calvis (:) + sol.sigman(:);

% The other way of calibrating the map, with A-team + diffuse emission present.
acm_cal = (sol.gainsol'*sol.gainsol) .* acc; 



antmask = zeros (size (acc));
posmask = zeros (size (poslocal));
rem_ants = length(acc) - length(flagant);
for ind = 1:length(flagant)
  antmask (flagant(ind), :) = 1; antmask (:,flagant(ind)) = 1;
  posmask (flagant(ind), :) = 1;
end
acc_fl = reshape (acm_cal(antmask ~= 1), [rem_ants, rem_ants]);
poslocal_fl = reshape(poslocal(posmask ~=1), [rem_ants, 3]);
posITRF_fl = reshape (posITRF(posmask ~= 1), [rem_ants, 3]);


% Show the UV plane density of the visibilities before and after weighting
gridnoweight = genvisden (posITRF_fl, 2, ones (1,prod(size(acc_fl)))); % No weighting
gridwithweight = genvisden (posITRF_FL, 2, 1./weight);
figure; mesh(gridnoweight); xlabel ('U(m)'); ylabel ('V(m)'); zlabel ('Visibility count (No weighting)');
figure; mesh(gridwithweight); xlabel ('U(m)'); ylabel ('V(m)'); zlabel ('Visibility count (With weighting)');

%% Create calibrated maps
innermap = acm2skyimage (acc_fl, poslocal_fl(:,1), poslocal_fl(:,2), fobs, l, l);
mask = NaN(size (innermap));
mask(meshgrid(l).^2 + meshgrid(l).'.^2 < 1) = 1;
imagesc (l, l, innermap.* mask);
colorbar;
set(gca, 'FontSize', 16);
title (sprintf ('LBA_INNER calib. Time: %s Freq: %.2f', datestr(mjdsec2datenum(tobs)), fobs));
axis equal
axis tight
set (gca, 'YDir', 'Normal'); % To match orientation with station images
set (gca, 'XDir', 'Reverse'); % To match orientation with station images

ylabel('South $\leftarrow$ m $\rightarrow$ North', 'interpreter', 'latex');
xlabel('East $\leftarrow$ l $\rightarrow$ West', 'interpreter', 'latex');
drawnow ();
if (doprint == 1)
    print (sprintf ('~/Documents/AARTFAAC_Project_SW_system_plan/afaac_psfshape/LBA_INNER_analysis/afaac_cal.png', stat), '-dpng', '-r300');
end;

% Create calibrated single station maps
for stat = 1:6
    stflag = setdiff ([1:288], setdiff ((stat-1)*48 + [1:48], flagant));
    antmask = zeros (size (acc));
    posmask = zeros (size (poslocal));
    rem_ants = length(acc) - length(stflag);
    for ind = 1:length(stflag)
        antmask (stflag(ind), :) = 1; antmask (:,stflag(ind)) = 1;
        posmask (stflag(ind), :) = 1;
    end
    acc_st = reshape (acm_cal(antmask ~= 1), [rem_ants, rem_ants]);
    poslocal_fl = reshape(poslocal(posmask ~=1), [rem_ants, 3]); 
    
    statmap = acm2skyimage (acc_st, poslocal_fl(:,1), poslocal_fl(:,2), fobs, l, l);
    imagesc (l, l, statmap.* mask);
    colorbar;
    set(gca, 'FontSize', 16);
    title (sprintf ('St. %d, calib Time: %s Freq: %.2f', stat, datestr(mjdsec2datenum(tobs)), fobs));
    axis equal
    axis tight
    set (gca, 'YDir', 'Normal'); % To match orientation with station images
    set (gca, 'XDir', 'Reverse'); % To match orientation with station images

    ylabel('South $\leftarrow$ m $\rightarrow$ North', 'interpreter', 'latex');
    xlabel('East $\leftarrow$ l $\rightarrow$ West', 'interpreter', 'latex');
    drawnow ();
    if (doprint == 1)
        print (sprintf ('~/Documents/AARTFAAC_Project_SW_system_plan/afaac_psfshape/LBA_INNER_analysis/st%d_cal.png', stat), '-dpng', '-r300');
    end;
end;


