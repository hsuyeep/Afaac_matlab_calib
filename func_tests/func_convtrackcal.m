% Script to generate solutions using a single convergent calibration, followed
% by tracking calibration using varying number of iterations. Needed to map 
% the convergence of tracking calibration to a convergent calibration solu.
% pep/24Feb13

addpath ../
fname = '~/WORK/AARTFAAC/Reobs/11Jul12/LBA_OUTER_BAND_SPREAD/SB004_LBA_OUTER_SPREAD_1ch.bin';

col = {'-bo', '-mo', '-ro', '-ko', '-go', '-yo', '-co',  ...
	   '-b*', '-m*', '-r*', '-k*', '-g*', '-y*', '-c*'};
flagant = [51, 238, 273]; 
visamphithresh= 1.5;% Reject visibilities with median >visampthresh*median.
visamplothresh= 0.5;% Reject visibilities with median >visampthresh*median.
debuglev = 0;
prevsol = [];
nrec = 120;           % Number of timeslices on which to operate.
nstefiter = 8;      % Maximum number of stefcal iterations in tracking cal.
ncalextiter = 4;    % Maximum number of cal_ext_stefcal iter. in tracking cal.
%err = zeros (nrec, nstefiter, ncalextiter);
err = zeros (nrec, nstefiter);
sigmas_track = zeros (nrec, 5);
sigmas_conv  = sigmas_track;

%{
% Carry out convergent calibration to generate initial solution for tracking.
prevsol = pelican_sunAteamsub (acc, tobs, freq, uvflag, ... 
			   	currflagant, debuglev, 1, [], []);

% Characteristics from convergent solution
fprintf (1, 'Conv.cal: calext_iter: %d, pinv_sol: %f, stefsol_iter: %d, stefsol.tol: %f\n', prevsol.calext_iters, prevsol.pinv_sol,...
prevsol.stefsol.iter,  prevsol.stefsol.tol);

[acc, tobs, freq] = readms2float (fid, -1, -1, 288);
[uvflag, missant] = flagdeadcorr (acc, tobs, freq, visamphithresh, ...
										visamplothresh);
%}

calplt = figure();
imgplt = figure();
imgparm = [];

for niter=1:1 %[1 2 4 8] % 1:nstefiter;
	fprintf (2, '\n\nTesting with stefsol iterations: %d\n', niter);
	clear pelican_tracking_cal;
	clear pelican_sunAteamsub;
	% err = zeros (1, nrec);
	% err = zeros (nrec, nstefiter);

	clf;
	fid = fopen (fname, 'rb'); % Reopen file for every test, for identical data.
	for ind = 1:50
		[acc, tobs, freq] = readms2float (fid, -1, -1, 288);
	end;
	% Update flagged visibilities
	[uvflag, missant] = flagdeadcorr (acc, tobs, freq, visamphithresh, ...
									visamplothresh);
	currflagant = unique ([flagant missant]);

	% Try out giving a conv. cal solu. as initial estimates for tracking.
	prevsol = pelican_sunAteamsub (acc, tobs, freq, uvflag, ... 
    				   	currflagant, debuglev, 1, [], []);
	% prevsol = [];
	for ts=1:nrec
  		fprintf (2, '\n\nTimeslice: %f\n', tobs);
		% Carry out convergent calib. to generate comparison vector.
		currconvsol = pelican_sunAteamsub (acc, tobs, freq, uvflag, ... 
							   	currflagant, debuglev, 1, [], []);
		sigmas_conv (ts, :) = currconvsol.sigmas;

		% Carry out tracking calib. for different number of stefcal iterations.
		currtracksol = pelican_tracking_cal (acc, tobs, freq, uvflag, ... 
							currflagant, debuglev, 1, prevsol, [], 8);
		% sigmas_track (ts, :) = currtracksol.sigmas;
	
		% Display this solution's characteristics.
		fprintf (1, ... 
		'Calext iters: %d, stefsol niters: %d, steftol: %f, pinv_sol: %f\n',...
				 currtracksol.calext_iters, currtracksol.stefsol.iter,...
				 currtracksol.stefsol.tol, currtracksol.pinv_sol);
	
		% Display estimated parameters for second timeslice from both calib.
		figure (calplt);
		subplot (2,2,1);
%		plot (real(currconvsol.gainsol) - real(currtracksol.gainsol), ...
%				char(col(mod(ts, length(col))+1)));
%		ylabel ('Gain real difference');
		plot (abs(currconvsol.gainsol) - abs (currtracksol.gainsol), ... 
			  char(col(mod(ts, length(col))+1))); 
		xlabel ('Ant. num');
		ylabel ('Gain ampl. difference');
		title (sprintf ('Time: %.2f, iter: %d, caliter: %d, stefiter: %d', ... 
		  tobs, niter, currtracksol.calext_iters, currtracksol.stefcal_iters));
		hold on;

		subplot (2,2,2);
		plot (angle(currconvsol.gainsol) - angle (currtracksol.gainsol), ... 
			  char(col(mod(ts, length(col))+1))); 
%		plot (imag(currconvsol.gainsol) - imag(currtracksol.gainsol), ...
%				char(col(mod(ts, length(col))+1)));
%		ylabel ('Gain imag difference');
		ylabel ('Gain phase difference');
		xlabel ('Ant. num');
		hold on;

		subplot (2,2,3);
		% Error norm, in dB
		% err (ts, niter) = 10*log10(norm (currconvsol.gainsol - ...
	 	% 						currtracksol.gainsol, 'fro'));
		% ylabel ('10*log10(||convgain - trackgain||Fro)');
		% title ('Residual Frobenius norms');

		% Relative error in percent
		err (ts, niter) = 100*sum(abs(currconvsol.gainsol - ...
	  	 			currtracksol.gainsol))/sum(abs(currconvsol.gainsol));
		% err (ts, niter) = 100*sum(angle(currconvsol.gainsol - ...
	  	% 			currtracksol.gainsol))/sum(angle(currconvsol.gainsol));
		plot (err (:, niter), '-ro');
		xlabel ('Timeslice number');
		ylabel ('Relative % error');
		title ('Error difference between convergent and tracking cal');

		subplot (2,2,4);
		plot (sigmas_conv (ts, :) - sigmas_track (ts, :), ...
			  char(col(mod(ts, length(col))+1))); 
		xlabel ('Source number');
		ylabel ('sigmas difference');
		leg {ts} = sprintf ('ts:%02d', ts);
		hold on;
		legend (leg);
		drawnow ();
		% pause;

		% --------- NOTE NOTE NOTE!!! This part decides the tracking! ---------
		% prevsol = currtracksol; % Always prev. tracking sols within timeslices.

		% Generate image statistics
%		calvis (:,:,1) = currconvsol.calvis;
%		calvis (:,:,2) = currtracksol.calvis;
%		[img, imgparmret] = ... 
%		  cmpimg (calvis, 2, imgparm, [currconvsol.tobs, currtracksol.tobs], ...
%				  [currconvsol.freq, currtracksol.freq], currflagant);
%		imgparm = imgparmret;
%		imgparm.dofft = 1;
%		figure(imgplt);
%		subplot (1,2,1);
%		plot ([1:ts], img.nlse(2), '-bo');
%		hold on;
%
%		subplot (1,2,2);
%		a = hist (reshape(img.diffmap(:,:,2), 1, 262144));
%		plot (a(2:end),'-ro');
%		hold on;

		% imagesc (img.diffmap(:,:,2));

		[acc, tobs, freq] = readms2float (fid, -1, -1, 288);
	end;
	fclose (fid);
	% pause;
end;

%{
% Plot errs 
figure;
for ind=0:size (err,1)-1
	% fprintf (1, '%d %s\n', ind, col{mod(ind, length(col)) + 1});
	plot (err (ind+1, :), char(col(mod(ind, length(col)) + 1)));
	leg {ind+1} = sprintf ('ts:%02d', ind + 1);
	hold on;
end;
legend (leg);
xlabel ('Number of iterations');
ylabel ('10*log10(||convgain - trackgain||Fro)');

figure;
for ind=0:size (err,2)-1
	plot (err (:, ind+1), char(col(mod(ind, length(col)) + 1)));
	leg {ind+1} = sprintf ('No. iter: %02d', ind+1);
	hold on;
end;
legend (leg);
xlabel ('Number of Timeslices');
ylabel ('10*log10(||convgain - trackgain||Fro)');
%}
