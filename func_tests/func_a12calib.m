% Script for trying out calibration of the A12 dataset.
% pep/24Nov15

% The problem: An A-12 dataset on regular calibration seems to calibrate, but not
% subtract out CasA correctly. This is because the flux of CasA is estimated
% as -ve by LS imaging. Note that the -ve LS imaging values themselves don't
% seem to be the cause since the first iteration of A6 calibration also results
% in -ve fluxes. The A-6 subset calibrates perfectly.

% Differences between A6 and A12:
% - W component of A12 ranges between 
% - Higher resolution. Max baseline is ~1Km, i.e., 20' at 60 MHz.
%   However, the A-team sources should be a max of a few arcmin, so should not
%   be a problem.
% - Non uniform UV coverage. Average visibility density is

% Approaches:
% 1 Calibrate A6 separately, the remaining 6 stations separately. 
%   Will probably not work due to issues of higher resolution remaining.
% 2 Calibrate A6 separately, apply the solutions to input ACM, then calibrate
%	the remaining stations.
% 3 Calibrate the remaining A-12 stations in subsets with A-6 stations, such
%   that the max baseline remains within 300 m.
% 4 Calibrate individual stations of the inner circle, apply calibration, then
% 	carry out a 12 station calibration.
% 5 Calibrate pairs of stations such that LS imaging output is not negative.


load ('/dop312_0/prasad/GPU_CORR_DAT/afaac-12/lbaouter_06Nv15/sb0_1447252930.vis_1447252960-1447252969.mat');
fobs = 195312.5*296; % Hz.
addpath '~/Documents/AARTFAAC_Project_SW_system_plan/afaac_GPU_interface/src'
flagant = [324, 373, 419, 492, 527, 537, 538];
[acm_t, tmjdsec,fobs,map,l] = gengpuimg (acm, 576,tobs,fobs,[1:15],[],[],[],0,0);

acm_uncal= zeros (1, 576, 576);
acm_uncal(1, :, :) = acm_t(1,:,:,1);

uvflag = eye(576);

%%%%% Approach 1.
tic; sola6  = pelican_sunAteamsub (conj(squeeze(acm_t(1,:,:,1))), tmjdsec(1), fobs, eye(576), union ([289:576], flagant), 0, 1,[], [], 'poslocal_afaac12_outer.mat', [], []);toc;
tic; sola12 = pelican_sunAteamsub (conj(squeeze(acm_t(1,:,:,1))), tmjdsec(1), fobs, eye(576), union ([  1:288], flagant), 0, 1, [], [], 'poslocal_afaac12_outer.mat', [], []); toc;

% Now plot the calibration solutions for the two sets of 6 stations.
figure;
subplot (2,1,1);
plot (abs(sola6.gainsol)); hold on; plot (abs(sola12.gainsol), 'r');
subplot (2,1,2);
plot (angle(sola6.gainsol)); hold on; plot (angle(sola12.gainsol), 'r');

% Generate calibrated images from the two sets of 6 stations (just to see how images from
% the outer ring look).buus
antmask = zeros (576);
flags = union ([289:576], flagant);
antmask(flags,:) = 1; antmask (:,flags) = 1;
acm_cal = zeros (1, 576, 576);
acm_tmp = zeros (576);
acm_tmp (antmask == 0) = sola6.calvis;
acm_cal (1, :, :) =  acm_tmp;
[l,m,calimga6,rdacc,locacc] = genfftimage(acm_cal, 1, 0,0,1,flags, 'poslocal_afaac12_outer.mat',0,[],0,0,0,0,0,tmjdsec(1), fobs);
imagesc (-l, -m, real(calimga6)); colorbar; title (sprintf ('A6-cal:%s', datestr(mjdsec2datenum(tmjdsec(1)))));

antmask = zeros (576);
flags = union ([  1:288], flagant);
antmask(flags,:) = 1; antmask (:,flags) = 1;
acm_cal = zeros (1, 576, 576);
acm_tmp = zeros (576);
acm_tmp (antmask == 0) = sola12.calvis;
acm_cal (1, :, :) =  acm_tmp;
[l,m,calimgcirc,rdacc,locacc] = genfftimage(acm_cal, 1, 0,0,1,flags, 'poslocal_afaac12_outer.mat',0,[],0,0,0,0,0,tmjdsec(1), fobs);
imagesc (-l, -m, real(calimgcirc)); colorbar; title (sprintf ('Inner circle:%s', datestr(mjdsec2datenum(tmjdsec(1)))));
% NOTE: Generated inner circle image is rubbish. Because the calibration did not
% go through.


% Apply the generated calibration vectors for all antennas to the data 
app1_gains = sola6.gainsol + sola12.gainsol;
acm_precal = (app1_gains'*app1_gains) .* conj (squeeze(acm_t(1,:,:,1)));

% Need to save solutions to disk to clear read only data in calibration routine
save ('/tmp/func_a12calib_approach1.mat', 'sola6', 'sola12', 'app1_gains', 'acm_precal', 'tmjdsec', 'fobs', 'flagant');
clear all; clear pelican_sunAteamsub;
load ('/tmp/func_a12calib_approach1.mat');


% Calibrate the ACM now precalibrated in subsets of 6 stations.
tic; solapp1 = pelican_sunAteamsub (acm_precal, tmjdsec(1), fobs, eye(576), flagant, 0, 1, [], [], 'poslocal_afaac12_outer.mat', [], []); toc;

% Show final A12 calibrated, Ateam subtracted image.
antmask = zeros (576);
flags = flagant;
antmask(flags,:) = 1; antmask (:,flags) = 1;
acm_cal = zeros (1, 576, 576);
acm_tmp = zeros (576);
acm_tmp (antmask == 0) = solapp1.calvis;
acm_cal (1, :, :) =  acm_tmp;
[l,m,calimga12,rdacc,locacc] = genfftimage(acm_cal, 1, 0,0,1,flags, 'poslocal_afaac12_outer.mat',0,[],0,0,0,0,0,tmjdsec(1), fobs);
imagesc (-l, -m, real(calimga12)); colorbar; title (sprintf ('Approach 1: A12 after precalibration:%s', datestr(mjdsec2datenum(tmjdsec(1)))));

%%%%%% Approach 2

% Apply the calibration vector of the 6 stations
sola6.gainsol (289:576) = 1;
acm_precal = (sola6.gainsol'*sola6.gainsol).*conj(squeeze(acm_t(1,:,:,1)));

% Now calibrate all 12 stations with the calibration vector applied.
tic; solapp2 = pelican_sunAteamsub (acm_precal, tmjdsec(1), fobs, eye(576), flagant, 0, 1, [], [], 'poslocal_afaac12_outer.mat', [], []); toc;

% Generate a calibrated, Ateam subtracted image.
antmask = zeros (576);
flags = flagant;
antmask(flags,:) = 1; antmask (:,flags) = 1;
acm_cal = zeros (1, 576, 576);
acm_tmp = zeros (576);
acm_tmp (antmask == 0) = solapp1.calvis;
acm_cal (1, :, :) =  acm_tmp;
[l,m,calimgapp2,rdacc,locacc] = genfftimage(acm_cal, 1, 0,0,1,flags, 'poslocal_afaac12_outer.mat',0,[],0,0,0,0,0,tmjdsec(1), fobs);
imagesc (-l, -m, real(calimgapp2)); colorbar; title (sprintf ('Approach 2: A12 after precalibration:%s', datestr(mjdsec2datenum(tmjdsec(1)))));
% NOTE: Same result as approach 1 after precalibration. Looks like the outer
% antennas are not calibrating at all.

%%%%%% Approach 3

% Determine subsets of stations within a fixed radius some station in the
% superterp
% Distances of each station from CS002, by taking the first dipole in the
% antennaset. (TODO: Check if this is the central one!)
uloc =  meshgrid (poslocal_outer(:,1)) - meshgrid (poslocal_outer(:,1).';
vloc =  meshgrid (poslocal_outer(:,2)) - meshgrid (poslocal_outer(:,2).';
uvdist = sqrt (uloc.^2 + vloc.^2);
station_dips = [0, 6:11]*48 + 1;
uvdist_stations = uvdist(station_dips); % wont work, is a matrix
dist = sqrt(poslocal_outer(station_dips,1).^2 + poslocal_outer(station_dips,2)^2);

% Find stations within 300m from each other
good_stations = find (uvdist_stations < 300);

% Calibrate by adding those stations to the superterp, if possible, or creating
% subsets of stations with max baselines within 300m.

%%%%%%%%%%% Approach 4
% Calibrate individual stations of the inner circle of stations first.
tic; solcs02= pelican_sunAteamsub (squeeze(acm_t(1,:,:,1)), tmjdsec(1), fobs, eye(576), union (flagant, setdiff([1:576], [1:48])), 0, 1, [], [], 'poslocal_afaac12_outer.mat', [], []); toc;
tic; solcs11= pelican_sunAteamsub (squeeze(acm_t(1,:,:,1)), tmjdsec(1), fobs, eye(576), union (flagant, setdiff([1:576], [289:336])), 0, 1, [], [], 'poslocal_afaac12_outer.mat', [], []); toc;
tic; solcs13= pelican_sunAteamsub (squeeze(acm_t(1,:,:,1)), tmjdsec(1), fobs, eye(576), union (flagant, setdiff([1:576], [337:384])), 0, 1, [], [], 'poslocal_afaac12_outer.mat', [], []); toc;
tic; solcs17= pelican_sunAteamsub (squeeze(acm_t(1,:,:,1)), tmjdsec(1), fobs, eye(576), union (flagant, setdiff([1:576], [385:432])), 0, 1, [], [], 'poslocal_afaac12_outer.mat', [], []); toc;

tic; solapp3= pelican_sunAteamsub (squeeze(acm_t(1,:,:1)), tmjdsec(1), fobs, eye(576), setdiff([1:576], [289:336]), 0, 1, [], [], 'poslocal_afaac12_outer.mat', [], []); toc;
solcs11 = pelican_sunAteamsub (squeeze(acm_t(1,:,:1)), tmjdsec(1), fobs, eye(576), setdiff ([1:576], [337:384]), 0, 1, [], [], 'poslocal_afaac12_outer.mat', [], []); 
solcs13 = pelican_sunAteamsub (squeeze(acm_t(1,:,:1)), tmjdsec(1), fobs, eye(576), setdiff ([1:576], [385:432]), 0, 1, [], [], 'poslocal_afaac12_outer.mat', [], []); 
solcs17 = pelican_sunAteamsub (squeeze(acm_t(1,:,:1)), tmjdsec(1), fobs, eye(576), setdiff ([1:576], [433:480]), 0, 1, [], [], 'poslocal_afaac12_outer.mat', [], []); 
solcs21 = pelican_sunAteamsub (squeeze(acm_t(1,:,:1)), tmjdsec(1), fobs, eye(576), setdiff ([1:576], [481:528]), 0, 1, [], [], 'poslocal_afaac12_outer.mat', [], []); 
solcs32 = pelican_sunAteamsub (squeeze(acm_t(1,:,:1)), tmjdsec(1), fobs, eye(576), setdiff ([1:576], [529:576]), 0, 1, [], [], 'poslocal_afaac12_outer.mat', [], []); 
