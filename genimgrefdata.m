% Script to generate images from a reference dataset, at the location:
% REFDATA_A6_LBA_OUTER_1419000369_SB296_1419000369-1419000378.mat
% Used for debugging AARTFAAC-12 data.
% Start of subband SB296 = 195312.5*296 - 195312.5/2 = 57714843.750 Hz.
% With 64 channels, channel i has a central frequency of 195312.5*296 - (195312.5/2) + (i+1)*3051.757812 Hz
% pep/23Nov15

% load ('/dop312_0/prasad/GPU_CORR_DAT/afaac-12/lbaouter_06Nv15/sb0_1446772267.vis_1446772297-1446772298.mat');
% flagant = [84, 97:144, 179, 250, 262, 469, 492, 493, 495, 496, 538, 551];
close all; clear all;
addpath '~/Documents/AARTFAAC_Project_SW_system_plan/afaac_GPU_interface/src/'
load ('~/WORK/AARTFAAC/Afaac_matlab_calib/REFDATA_A6_LBA_OUTER_1419000369_SB296_1419000369-1419000378.mat');

nant = 288;
freq = 195312.5*296 - 195312.5/2;
flagant = [140,149,260];
[acm_t, tmjdsec,fobs,map,l] = gengpuimg (acm, nant,tobs,freq,[1:63],[],[],[],0,0);

acm_uncal= zeros (1, nant, nant);
acm_uncal(1, :, :) = acm_t(1,:,:,1);

% Uncalibrated image:
[l,m,uncalimg,rdacc,locacc] = genfftimage(conj(acm_uncal), 1, 0,0,1,flagant, 'poslocal_outer.mat',0,[],0,0,0,0,0,tmjdsec(1), freq);
imagesc (-l, -m, real(uncalimg)); colorbar; title (sprintf ('A6-uncal:%s', datestr(mjdsec2datenum(tmjdsec(1)))));
axis equal
axis tight
set (gca, 'YDir', 'Normal'); % To match orientation with station images
set (gca, 'XDir', 'Reverse'); % To match orientation with station images
ylabel('South \leftarrow m \rightarrow North'); % , 'interpreter', 'latex');
xlabel('East \leftarrow l \rightarrow West'); % , 'interpreter', 'latex');


% Calibrate A6 dipoles alone
tic; sola6 = pelican_sunAteamsub(conj(squeeze(acm_t(1,:, :,1))), tmjdsec(1), fobs, eye(nant), flagant, 0, 1,[], [], 'poslocal_outer.mat', [], []);toc;
antmask = zeros (nant);
antmask (flagant,:) = 1; antmask(:,flagant) = 1;
acm_cal = zeros (1, nant, nant);
acm_tmp = zeros (nant);
acm_tmp (antmask == 0) = sola6.calvis; 
acm_cal (1, :, :) =  acm_tmp+ eye(288); % To make all pixels positive.;
[l,m,calimga6,rdacc,locacc] = genfftimage(acm_cal, 1, 0,0,1, flagant, 'poslocal_outer.mat',0,[],0,0,0,0,0,tmjdsec(1), freq);
imagesc (real(calimga6)); title ('Calibrated A6');

% Generate images without A-team subtraction
acm_cal (1,:,:) = (sola6.gainsol'*sola6.gainsol) .* squeeze(acm_t(1,:,:,1));
[l,m,calimga6_ateam,rdacc,locacc] = genfftimage(acm_cal, 1, 0,0,1, flagant, 'poslocal_outer.mat',0,[],0,0,0,0,0,tmjdsec(1), freq);
imagesc (real(calimga6_ateam)); title ('Calibrated A6, with Ateam');
%save (sprintf ('a6stuff_%s.mat', datestr (now, 30)));
