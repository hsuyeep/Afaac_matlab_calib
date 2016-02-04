% Script to generate images from a .mat file, as generated by gpu2mat.py
% Used for debugging AARTFAAC-12 data.
% pep/09Nov15

% load ('/dop312_0/prasad/GPU_CORR_DAT/afaac-12/lbaouter_06Nv15/sb0_1446772267.vis_1446772297-1446772298.mat');
% flagant = [84, 97:144, 179, 250, 262, 469, 492, 493, 495, 496, 538, 551];
close all; clear all;
addpath '~/Documents/AARTFAAC_Project_SW_system_plan/afaac_GPU_interface/src/'

% Calibrate a 'good; a-6 observation.
% fid = fopen ('~/WORK/AARTFAAC/Reobs/20Nov13/r02/SB002_LBA_OUTER_8b2sbr02_1ch.bin', 'rb');
fid = fopen ('~/WORK/AARTFAAC/Reobs/20Nov13/r03/SB002_LBA_OUTER_8b2sbr03_1ch_1337-1850.bin', 'rb');
[oldacm, oldtobs,oldfreq] = readms2float (fid, 10, -1, 288);
flagant = [129, 140, 149];
oldacm_tmp = zeros (1, 288, 288);
oldacm_tmp (1,:,:) = oldacm./sqrt((diag(oldacm)*diag(oldacm).')); % Whitened
ACM.
[l,m,uncaloldimg,rdacc,locacc] = genfftimage(oldacm_tmp, 1, 0,0,1, [], 'poslocal_outer.mat',0,[],0,0,0,0,0,oldtobs, oldfreq);
[dr_uncal, sig_uncal] = getimagedr (uncaloldimg, [], []);
fprintf (2, '<-- UnCalibrated image A6 (%s): SNR: %f, variance: %f\n', datestr(mjdsec2datenum(oldtobs)), dr_uncal, sig_uncal);

oldsol = pelican_sunAteamsub(oldacm, oldtobs, oldfreq, eye(288), flagant, 0, 1,[], [], 'poslocal_outer.mat', [], []);
acc =zeros(288);
mask = acc;
mask(flagant,:) = 1; mask(:,flagant) = 1;
acc(mask == 0) = oldsol.calvis;
oldacm_tmp (1,:,:) =acc;
[l,m,caloldimg,rdacc,locacc] = genfftimage(oldacm_tmp, 1, 0,0,1, [], 'poslocal_outer.mat',0,[],0,0,0,0,0,oldtobs, oldfreq);
[dr_cal, sig_cal] = getimagedr (caloldimg, [], []);
fprintf (2, '<-- Calibrated image A6 (%s): SNR: %f, variance: %f\n', datestr(mjdsec2datenum(oldtobs)), dr_cal, sig_cal);





% load ('/dop312_0/prasad/GPU_CORR_DAT/afaac-12/lbaouter_06Nv15/sb0_1447251401.vis_1447251432-1447251441.mat');
load ('/dop312_0/prasad/GPU_CORR_DAT/afaac-12/lbaouter_06Nv15/sb0_1447252930.vis_1447252960-1447252969.mat');
flagant = [324, 373, 419, 492, 527, 537, 538];

% Multiband observations (NOTE: Two stations have very low power)
% load ('/dop312_0/prasad/GPU_CORR_DAT/afaac-12/lbaouter_20Nov15_multiband/sb0_1448009421.vis_1448009451-1448009460.mat');
% flagant = [18, 84, 168, 262, 324, 373, 419, 492, 527, 537, 538];
[acm_t, tmjdsec,fobs,map,l] = gengpuimg (acm, 576,tobs,57217407.226562,[1:15],[],[],[],0,0);

acm_uncal= zeros (1, 576, 576);
tmp = squeeze(acm_t(1,:,:,1));
acm_uncal(1, :, :) = tmp./sqrt(diag(tmp)*diag(tmp).');

% Uncalibrated image: A12
[l,m,uncalimg,rdacc,locacc] = genfftimage(acm_uncal, 1, 0,0,1,flagant, 'poslocal_afaac12_outer.mat',0,[],0,0,0,0,0,tmjdsec(1), fobs);
figure;
imagesc (-l, -m, real(uncalimg)); colorbar; title (sprintf ('A12-uncal:%s', datestr(mjdsec2datenum(tmjdsec(1)))));
axis equal
axis tight
set (gca, 'YDir', 'Normal'); % To match orientation with station images
set (gca, 'XDir', 'Normal'); % To match orientation with station images
ylabel('South \leftarrow m \rightarrow North'); % , 'interpreter', 'latex');
xlabel('East \leftarrow l \rightarrow West'); % , 'interpreter', 'latex');
[dr_uncala12, sig_uncala12] = getimagedr (uncalimg, [], []);
fprintf (2, '<-- UnCalibrated image A12 (%s): SNR: %f, variance: %f\n', datestr(mjdsec2datenum(tmjdsec(1))), dr_uncala12, sig_uncala12);
text (-1, 0.9, sprintf ('SNR: %.2f, std: %.2f', dr_uncala12, sig_uncala12), 'Color', [1,1,1], 'Interpreter', 'latex', 'Fontsize', 16);
saveas (gcf, 'A12_uncal.png', 'png');

% Uncalibrated image: A6
tmp = squeeze(acm_t(1,1:288,1:288,1));
acm_tmp = zeros (1,288,288);
acm_tmp(1,:,:) = tmp./sqrt(diag(tmp)*diag(tmp).');
[l,m,uncalimga6,rdacc,locacc] = genfftimage(acm_tmp, 1, 0,0,1, [], 'poslocal_outer.mat',0,[],0,0,0,0,0,tmjdsec(1), fobs);
figure;
imagesc (-l, -m, real(uncalimga6)); colorbar; title (sprintf ('A6-uncal:%s', datestr(mjdsec2datenum(tmjdsec(1)))));
axis equal;
axis tight;
set (gca, 'YDir', 'Normal'); % To match orientation with station images
set (gca, 'XDir', 'Normal'); % To match orientation with station images
ylabel('South \leftarrow m \rightarrow North'); % , 'interpreter', 'latex');
xlabel('East \leftarrow l \rightarrow West'); % , 'interpreter', 'latex');
[dr_uncala6, sig_uncala6] = getimagedr (uncalimga6, [], []);
fprintf (2, '<-- UnCalibrated image A6 (%s): SNR: %f, variance: %f\n', datestr(mjdsec2datenum(tmjdsec(1))), dr_uncala6, sig_uncala6);
text (-1, 0.9, sprintf ('SNR: %.2f, std: %.2f', dr_uncala6, sig_uncala6), 'Color', [1,1,1], 'Interpreter', 'latex');
saveas (gcf, 'A6_uncal.png', 'png');

% Calibrate A6 dipoles alone
% sola6 = pelican_sunAteamsub(conj(squeeze(acm_t(1,1:288,1:288,1))), tmjdsec(1), fobs, eye(288), [], 0, 1,[], [], 'poslocal_outer.mat', [], []);
flagant_288 = [289:576];
tic; sola6 = pelican_sunAteamsub(conj(squeeze(acm_t(1,:, :,1))), tmjdsec(1), fobs, eye(576), flagant_288, 0, 1,[], [], 'poslocal_afaac12_outer.mat', [], []);toc;
antmask = zeros (576);
antmask (flagant_288,:) = 1; antmask(:,flagant_288) = 1;
acm_cal = zeros (1, 576, 576);
acm_tmp = zeros (576);
acm_tmp (antmask == 0) = sola6.calvis;
acm_cal (1, :, :) =  acm_tmp;
[l,m,calimga6,rdacc,locacc] = genfftimage(acm_cal, 1, 0,0,1, flagant_288, 'poslocal_afaac12_outer.mat',0,[],0,0,0,0,0,tmjdsec(1), fobs);
figure();
imagesc (real(calimga6)); colorbar(); title ('Calibrated A6');
[dr_cala6, sig_cala6] = getimagedr (calimga6, [], []);
fprintf (2, '<-- Calibrated image A6 (%s): SNR: %f, variance: %f\n', datestr(mjdsec2datenum(mjdsec(1))), dr_cala6, sig_uncala6);

% Generate images without A-team subtraction
acm_cal (1,:,:) = (sola6.gainsol'*sola6.gainsol) .* squeeze(acm_t(1,:,:,1));
[l,m,calimga6_ateam,rdacc,locacc] = genfftimage(acm_cal, 1, 0,0,1, flagant_288, 'poslocal_afaac12_outer.mat',0,[],0,0,0,0,0,tmjdsec(1), 60000000);
figure;
imagesc (real(calimga6_ateam)); colorbar(); title ('Calibrated A6, with Ateam');
[dr_cala6ateam, sig_cala6ateam] = getimagedr (calimga6_ateam, [], []);
fprintf (2, '<-- Calibrated image A6+Ateam(%s): SNR: %f, variance: %f\n', datestr(mjdsec2datenum(tmjdsec(1))), dr_cala6ateam, sig_cala6ateam);
%save (sprintf ('a6stuff_%s.mat', datestr (now, 30)));


%%%%%%%%%%% Calibrate A12!
tic; 
sol = pelican_sunAteamsub(conj(squeeze(acm_t(1,:,:,1))), tmjdsec(1), fobs, eye(576), flagant, 0, 1,[], [], 'poslocal_afaac12_outer.mat', [], []); 
toc;

%sol_6stat = pelican_sunAteamsub(conj(squeeze(acm_t(1,:,:,1))), tmjdsec(1), fobs, eye(576), union(flagant, [289:576]), 0, 1,[], [], 'poslocal_afaac12_outer.mat', [], []); 
%sol_7stat = pelican_sunAteamsub(conj(squeeze(acm_t(1,:,:,1))), tmjdsec(1), fobs, eye(576), union(flagant, [337:576]), 0, 1,[], [], 'poslocal_afaac12_outer.mat', [], []); 
%sol_8stat = pelican_sunAteamsub(conj(squeeze(acm_t(1,:,:,1))), tmjdsec(1), fobs, eye(576), union(flagant, [385:576]), 0, 1,[], [], 'poslocal_afaac12_outer.mat', [], []); 
%sol_9stat = pelican_sunAteamsub(conj(squeeze(acm_t(1,:,:,1))), tmjdsec(1), fobs, eye(576), union(flagant, [433:576]), 0, 1,[], [], 'poslocal_afaac12_outer.mat', [], []); 
%sol_10stat = pelican_sunAteamsub(conj(squeeze(acm_t(1,:,:1))), tmjdsec(1), fobs, eye(576), union(flagant, [481:576]), 0, 1,[], [], 'poslocal_afaac12_outer.mat', [], []); 
%sol_11stat = pelican_sunAteamsub(conj(squeeze(acm_t(1,:,:1))), tmjdsec(1), fobs, eye(576), union(flagant, [529:576]), 0, 1,[], [], 'poslocal_afaac12_outer.mat', [], []); 
antmask = zeros (576);
antmask(flagant,:) = 1; antmask (:,flagant) = 1;
acm_cal = zeros (1, 576, 576);
acm_tmp = zeros (576);
acm_tmp (antmask == 0) = sol.calvis;
acm_cal (1, :, :) =  acm_tmp;
[l,m,calimg,rdacc,locacc] = genfftimage(acm_cal, 1, 0,0,1,flagant, 'poslocal_afaac12_outer.mat',0,[],0,0,0,0,0,tmjdsec(1), 60000000);
figure;
imagesc (-l, -m, real(calimg)); colorbar();  title (sprintf ('A12-cal with Ateam:%s', datestr(mjdsec2datenum(tmjdsec(1)))));
axis equal
axis tight
set (gca, 'YDir', 'Normal'); % To match orientation with station images
set (gca, 'XDir', 'Reverse'); % To match orientation with station images
ylabel('South \leftarrow m \rightarrow North'); % , 'interpreter', 'latex');
xlabel('East \leftarrow l \rightarrow West'); % , 'interpreter', 'latex');

% Generate images without A-team subtraction
acm_cal (1,:,:) = (sol.gainsol'*sol.gainsol) .* squeeze(acm_t(1,:,:,1));
[l,m,calimg_ateam,rdacc,locacc] = genfftimage(acm_cal, 1, 0,0,1, flagant, 'poslocal_afaac12_outer.mat',0,[],0,0,0,0,0,tmjdsec(1), 60000000);
figure;
imagesc (-l, -m, real(calimg_ateam)); colorbar(); title (sprintf ('A12-cal with Ateam:%s', datestr(mjdsec2datenum(tmjdsec(1)))));


%uncalhdl = figure();
%imagesc (-l, -m, real(uncalimg)); colorbar; title (sprintf ('A12-uncalibrated image:%s', tmjdsec));
%title (sprintf ('Time: %s, A12-uncal', datestr(mjdsec2datenum(tmjdsec(1)))));
%axis equal
%axis tight
%% set (gca, 'YDir', 'Normal'); % To match orientation with station images
%% set (gca, 'XDir', 'Reverse'); % To match orientation with station images
%
%ylabel('South $\leftarrow$ m $\rightarrow$ North', 'interpreter', 'latex');
%xlabel('East $\leftarrow$ l $\rightarrow$ West', 'interpreter', 'latex');
%load 'srclist3CR.mat';
%[sel, sel_l, sel_m] = overplotcat (tmjdsec(1), srclist3CR, 50, uncalhdl, 1);
%
%calhdl = figure();
%imagesc (-l, -m, real(calimg)); colorbar; title (sprintf ('A12-calibrated image:%s', tmjdsec));
%title (sprintf ('Time: %s, A12-cal', datestr(mjdsec2datenum(tmjdsec(1)))));
%axis equal
%axis tight
%[sel, sel_l, sel_m] = overplotcat (tmjdsec(1), srclist3CR, 50, uncalhdl, 1);
%
%subplot (212);
%imagesc (real(calimga6)); colorbar;title (sprintf ('A6-calibrated image:%s', tmjdsec));
