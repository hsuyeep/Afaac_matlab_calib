% Script to examine the effect of real data on tracking calibration in order to
% determine the effect of the number of iterations and the initial estimate of 
% the solutions on the final solution.
% pep/01Feb13

function [sol] = testtrackcal ()
	addpath '../';
	fid = fopen ('~/WORK/AARTFAAC/Reobs/11Jul12/LBA_OUTER_BAND_SPREAD/SB002_LBA_OUTER_SPREAD_1ch.bin', 'rb');

	% For LBA_OUTER_BAND_SPREAD, 18min data
    flagant = [51, 238, 273]; 
	debuglev = 2;
	ptSun = 1;
	visamphithresh= 1.5;% Reject visibilities with median >visampthresh*median.
	visamplothresh= 0.5;% Reject visibilities with median >visampthresh*median.
	
	ntimes = 10;
	[ac, t_obs, freq] = readms2float (fid, -1, -1, 288);

	[uvflag, missant] = flagdeadcorr (ac, t_obs, freq, visamphithresh, ...
										visamplothresh);
	flagant = unique ([flagant missant]);
	
	acc = zeros (size (ac, 1), size (ac, 2), ntimes);
	tobs = zeros (ntimes);
	freq = zeros (ntimes);		

	% Fill in raw data
	for ind = 1:ntimes
		[acc(:,:,ind), tobs(ind), freq(ind)] = readms2float (fid, -1, -1, 288);
	end;

	ind = 1;
	% calibrate to convergence to generate truth value.
	convsol = pelican_sunAteamsub (acc (:,:,ind), t_obs(ind), freq(ind), ... 
							uvflag, flagant, debuglev, ptSun, 1, 30);
	clear pelican_sunAteamsub;
	clear first_call;

	% calibrate with tracking .
	tracksol = pelican_tracking_cal (acc (:,:,ind), t_obs(ind), freq(ind), ... 
							uvflag, flagant, debuglev, ptSun, 1, 30);
	% Now plot estimated parameters.
	% Cal_ext related, for this timeslice:
	
