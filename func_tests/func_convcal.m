% Script to test convergent calibration, while extracting out various solution
% characteristics.

addpath ../
fname = '~/WORK/AARTFAAC/Reobs/11Jul12/LBA_OUTER_BAND_SPREAD/SB004_LBA_OUTER_SPREAD_1ch.bin';

flagant = [51, 238, 273]; 
visamphithresh= 1.5;% Reject visibilities with median >visampthresh*median.
visamplothresh= 0.5;% Reject visibilities with median >visampthresh*median.
debuglev = 1;
prevsol = [];
nrec = 40;           % Number of timeslices on which to operate.

fid = fopen (fname, 'rb'); % Reopen file for every test, for identical data.
[acc, tobs, freq] = readms2float (fid, -1, -1, 288);

% Update flagged visibilities
[uvflag, missant] = flagdeadcorr (acc, tobs, freq, visamphithresh, ...
								visamplothresh);
currflagant = unique ([flagant missant]);

% Try out giving a conv. cal solu. as initial estimates for tracking.
currsol = pelican_sunAteamsub (acc, tobs, freq, uvflag, ... 
					   	currflagant, debuglev, 1, [], []);


