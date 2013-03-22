% Script to test the correctness of WSF wrt.
% 1. Simulated visibilities with varying amount of noise
% 2. Simulated visibilities with varying number of fminsearch iterations.
% 3. Real data with varying number of fminsearch iterations.
% 4. Long term real data with bad timeslices.
% pep/19Feb13

function [sol] = func_wsfsrcpos ()
	addpath ('../');

	% ----- Initialize RO data -----
	rodata.C       = 299792458;         % speed of light, m/s
	rodata.lon     = 6.869837540;       % longitude of CS002 in degrees
	rodata.lat     = 52.915122495;      % latitude of CS002 in degrees
	rodata.Nelem   = 288;               % Max. number of elements
	rodata.normal  = [0.598753, 0.072099, 0.797682].'; % Normal to CS002
	disp ('Loading 3CR catalog and local antenna positions.');
	load ('poslocal.mat', 'posITRF', 'poslocal'); 
	rodata.posITRF = posITRF;
	rodata.poslocal= poslocal;
	load srclist3CR;
	rodata.catalog = srclist3CR;  		% Sky catalog to use.

	% A-team from 3CR catalog for defining the model sky
	% Nsrc x 1 vector with source indices. -->NOTE<--: use 0 for the Sun
	rodata.srcsel =  [324, 283, 88, 179, 0];  % NOTE: No Sun being used!!


	% ----- Initialize calibration and imaging parameters ----- 
    calim.restriction    = 10; % Avoid vis. below 'restriction' 
                               % wavelengths length.
    calim.maxrestriction = 60; % Avoid vis. above 'maxrestriction' 
                               % meters (NOT wavelengths!)
    calim.debug = 0;      	   % Set debug level.
	calim.opt = optimset();    % Optimization parameters.
    calim.rem_ants = rodata.Nelem;
    calim.flagant = [];        % NOTE: Specify flagged antennas here!
	calim.uvflag = eye(rodata.Nelem);

	uloc = meshgrid (rodata.poslocal(:,1)) - ... 
			meshgrid (rodata.poslocal (:,1)).';
	vloc = meshgrid (rodata.poslocal(:,2)) - ... 
			meshgrid (rodata.poslocal (:,2)).';

	Nuv = 500; %1000        % size of gridded visibility matrix
	uvpad = 512; %1024      % specifies if any padding needs to be added

    % Calibration stopping conditions
    calim.diffstop = 1e-3;         % difference bet. calib. solutions
    calim.maxiter  = 10;            % Max. cal_ext iterations.
    calim.maxiter_gainsolv = 30;    

	% Flagging related
    antmask = zeros (rodata.Nelem);
    posmask = zeros (size (rodata.posITRF));
    if (length(calim.flagant) ~= 0)   
      	disp (['Flagant: ' num2str(calim.flagant)]);
    	calim.rem_ants = rodata.Nelem - length(calim.flagant);
    	for ind = 1:length(calim.flagant)
    		antmask (calim.flagant(ind), :) = 1; 
			antmask (:,calim.flagant(ind)) = 1;
    		posmask (calim.flagant(ind), :) = 1;
    	end
    	rodata.posITRF_fl = reshape(rodata.posITRF(posmask ~=1), ... 
        							[calim.rem_ants, 3]);
    	calim.uvflag  = reshape (calim.uvflag(antmask ~= 1), ... 
        				   [calim.rem_ants, calim.rem_ants]);

		% ---- Generate flagged uv coordinates. ----
		u=meshgrid(rodata.posITRF_fl(:, 1))-meshgrid(rodata.posITRF_fl(:, 1)).';
	   	v=meshgrid(rodata.posITRF_fl(:, 2))-meshgrid(rodata.posITRF_fl(:, 2)).';
	    w=meshgrid(rodata.posITRF_fl(:, 3))-meshgrid(rodata.posITRF_fl(:, 3)).';
	    uvw = [u(:), v(:), w(:)];
	    calim.uvdist = sqrt(sum(uvw.^2, 2) - (uvw * rodata.normal).^2);
    else
      	rodata.posITRF_fl = rodata.posITRF;
		% ---- Generate flagged uv coordinates. ----
		u=meshgrid(rodata.posITRF_fl(:, 1))-meshgrid(rodata.posITRF_fl(:, 1)).';
	   	v=meshgrid(rodata.posITRF_fl(:, 2))-meshgrid(rodata.posITRF_fl(:, 2)).';
	    w=meshgrid(rodata.posITRF_fl(:, 3))-meshgrid(rodata.posITRF_fl(:, 3)).';
	    uvw = [u(:), v(:), w(:)];
	    calim.uvdist = sqrt(sum(uvw.^2, 2) - (uvw * rodata.normal).^2);
    end
	

	col = {'-bo', '-mo', '-ro', '-ko', '-go', '-yo', '-co',... 
		   '-b*', '-m*', '-r*', '-k*', '-g*', '-y*', '-c*'};
	% Storage for estimated parameters
	simul = 0;
	if (simul == 1) fprintf (2, 'Simulation on\n'); end;
	maxiter = 20;
	thsrc_wsf  = zeros (length (rodata.srcsel), maxiter);
	phisrc_wsf = thsrc_wsf;
	thsrc_cat  = thsrc_wsf;
	phisrc_cat = thsrc_wsf;
	wsf_output = zeros (3, maxiter);
	iters = zeros (1, maxiter);

	% ----- Input parameters for simulation ----- 
	if (simul == 1)
		tobs = 4848774150.503317; % MJD sec, taken from LBA_OUTER_BAND_SPREAD obs.
		t_int  = 1; % 1 sec. integration
		f_int = 90000; % 90 KHz integration
		tobs_jd = tobs/86400.+2400000.5; % Convert to JD.
		freq = 69912719.726562;  % Hz
		duv = rodata.C/(freq*2); % Always Nyquist sample.
		abovejy = 300; 			 % Flux cutoff for generating a simulated sky.
	end;
	flagant = calim.flagant; % Flagged antennas.

	% Control the instantaneous output voltage from antennas. CURRENTLY UNUSED!
	antvolt = complex(ones(1, rodata.Nelem), zeros (1, rodata.Nelem));

	% This defines the variance of the gaussian noise to be added to the 
	% simulated gains per antenna. Currently identical, and random Tsys per
	% antenna.
	% antgainsig  = 0.1*ones (1, rodata.Nelem);
	antgainsig  = zeros (1, rodata.Nelem);

	% Control gains of the antennas, defined as a random complex number per
	% antenna on top of a constant station gain. Mimmics real situation with
	% gains being more or less calibrated within a station. These are the 
	% parameters being solved for by the calibration.

	% 1. No antenna gains at all.
	% fprintf (1, 'No antenna gains applied.\n');
	% stationgains = [0*exp(0*i*0.3), 0*exp(0*i*2.4), 0*exp(0*i*-0.0),...
 	% 			    0*exp(0*i*-0.2), 0*exp(0*i*2.8), 0*exp(0*i*-0.8)];

	% 2. Diverse antenna gains
%	 fprintf (1, 'Diverse antenna gains applied, diff. amp/ph per station.\n');
%	stationgains = [1.2*exp(1*i*1.3), 0.8*exp(1*i*2.4), 0.5*exp(1*i*-1.1),...
%  				    1.5*exp(1*i*-1.2), 1.8*exp(1*i*2.8), 0.7*exp(1*i*-1.8)];

	% 3. Only amplitude in antenna gains, no phases.
	% fprintf (1, 'Only ampl. in antenna gains applied, no phases.\n');
	%stationgains = [1.2*exp(0*i*1.3), 0.8*exp(0*i*2.4), 0.5*exp(0*i*-1.1),...
  	%			    1.5*exp(0*i*-1.2), 1.8*exp(0*i*2.8), 0.7*exp(0*i*-1.8)];
	
	% 4. Constant 1 as amp. gain, no phases.
	 fprintf (1, 'Const. 1 as ampl. gains applied, no phases.\n');
	 stationgains = [1*exp(0*i*1.3), 1*exp(0*i*2.4), 1*exp(0*i*-1.1),...
	 			    1*exp(0*i*-1.2), 1*exp(0*i*2.8), 1*exp(0*i*-1.8)];

	% 5. Constant 1 as amp. gain with diverse phases
	%fprintf (1, 'Const. 1 as ampl. gains applied with diverse phases.\n');
	% stationgains = [1*exp(1*i*1.3), 1*exp(1*i*2.4), 1*exp(1*i*-1.1),...
 	%			    1*exp(1*i*-1.2), 1*exp(1*i*2.8), 1*exp(1*i*-1.8)];

	
	antgains = complex (zeros (1, rodata.Nelem), zeros (1, rodata.Nelem));
	for ind = 1:length (stationgains)
		% Add a random component on station gain to get actual antenna gain.
		antset = [(ind-1)*48+(1:48)];
		re = [real(stationgains(ind)) + (antgainsig(antset).* randn (1,48))];
		im = [imag(stationgains(ind)) + (antgainsig(antset).* randn (1,48))];
		antgains (antset) = complex (re, im);
	end;

	% 1. All elements have identical, but diagonal noise.
	% sigman = 1.2*eye (rodata.Nelem); 
	sigman = eye (rodata.Nelem); 

	% 2. All elements have different, but diagonal noise.
	% sigman = diag(rand(1, rodata.Nelem));
	
	% 3. Visibilities < uvdist have a non-diagonal contribution. 
	% NOTE: This simulates only mutual coupling between antennas making up the
	% interferometer, and not diffuse emission.
%{
	uvdistlimit = 10; % NOTE: In lambda@freq. units
	mask = (reshape (calim.uvdist, [calim.rem_ants, calim.rem_ants]) < ... 
					   uvdistlimit * (rodata.C / freq));
	sigman = zeros (rodata.Nelem);
	for iind = 1:size (mask, 1)
		sigman (iind, iind) = randn (1, 1); % Autocorrelations
		for jind = iind+1:size (mask, 2) % Upper triangle only
			if (mask(iind, jind) == 1)
				sigman (iind, jind) = randn(1,1) + 1*i*randn(1,1);
				sigman (jind, iind) = conj (sigman (iind, jind));
			end;
		end;
	end;
%}

	% 4. Extreme case: Elements have different, and non-diagonal noise.
	% sigman = rand (rodata.Nelem);


	if (simul == 1)
		% ----- Generation of simulated sky -----

   	 [simsky_acc, simsky_A, mod_acc, mod_A, mod_up] = ... 
			genmodelsky (tobs, freq, [], rodata, ...
				calim, sigman, antgains, abovejy, calim.debug);
	 % cal1 = 1 ./ antgains;
	else
		% ----- Uncomment if input data from actual obs is needed -----
		% For LBA_OUTER_BAND_SPREAD, 18min data
    	flagant = [51, 238, 273]; 
 		fname = '~/WORK/AARTFAAC/Reobs/11Jul12/LBA_OUTER_BAND_SPREAD/SB004_LBA_OUTER_SPREAD_1ch.bin';
		fid = fopen (fname, 'rb');
		for ind = 1:12
			[simsky_acc, tobs, freq] = readms2float (fid, -1, -1, 288);
		end;
		fprintf (1, 'Time: %.2f, freq: %.2f\n', tobs, freq);

 		% ---- In case of flagged antennas, generate reshaped ACM ---- 
	    if (length(flagant) ~= 0)   
	      	disp (['Flagant: ' num2str(flagant)]);
	    	antmask = zeros (size (simsky_acc));
	    	posmask = zeros (size (rodata.posITRF));
	    	calim.rem_ants = length(simsky_acc) - length(flagant);
	    	for ind = 1:length(flagant)
	    		antmask (flagant(ind), :) = 1; antmask (:,flagant(ind)) = 1;
	    		posmask (flagant(ind), :) = 1;
	    	end
	    	simsky_acc = reshape (simsky_acc(antmask ~= 1), [calim.rem_ants, calim.rem_ants]);
	    	rodata.posITRF_fl = reshape(rodata.posITRF(posmask ~=1), ... 
	        							[calim.rem_ants, 3]);
	    	calim.uvflag  = reshape (calim.uvflag(antmask ~= 1), ... 
	        				   [calim.rem_ants, calim.rem_ants]);
	    	disp (['NOTE: ACM resized after flagging to ' num2str(size (simsky_acc))]);
			u=meshgrid(rodata.posITRF_fl(:, 1))-meshgrid(rodata.posITRF_fl(:, 1)).';
	   		v=meshgrid(rodata.posITRF_fl(:, 2))-meshgrid(rodata.posITRF_fl(:, 2)).';
	    	w=meshgrid(rodata.posITRF_fl(:, 3))-meshgrid(rodata.posITRF_fl(:, 3)).';
	    	uvw = [u(:), v(:), w(:)];
	    	calim.uvdist = sqrt(sum(uvw.^2, 2) - (uvw * rodata.normal).^2);
			% disp (['Rank/rcond of ACM: ' num2str(rank(acc)) ' ' num2str(rcond(acc))]);	
		end;
		tobs_jd = tobs/86400. + 2400000.5;
%   		[cal1, sigmas1, sigman] = statcal_stefcal (simsky_acc, tobs_jd, ... 
%									freq, rodata, calim, calim.uvflag);
	end;
    % whitening of the array covariance matrix (required for DOA estimation)
    % tic
    simsky_acc = simsky_acc ./ sqrt(diag(simsky_acc) * diag(simsky_acc).');

	% Carry out initial calibration of data
   	[cal1, sigmas1, sigman] = statcal_stefcal (simsky_acc, tobs_jd, ... 
									freq, rodata, calim, calim.uvflag);

    disp('First round of antenna calibration completed, sigmas1:'); 
    disp (sigmas1');

	% Generate initial estimates of positions of model sources
    if (sum(rodata.srcsel == 0) ~= 0) % Include Solar model
        [raSun, decSun] = SunRaDec(tobs_jd);                          % J2000
        rasrc = [rodata.catalog(rodata.srcsel(rodata.srcsel ~= 0)).alpha, ...
    			 raSun].';   % B1950
        decsrc = [rodata.catalog(rodata.srcsel(rodata.srcsel ~= 0)).delta, ...
    			  decSun].'; % B1950
        epoch = [true(length(rodata.srcsel(rodata.srcsel ~= 0)), 1); false];
        % disp (['Sun included, len of epoch ' num2str(length(epoch))]);
    else
        rasrc = [rodata.catalog(rodata.srcsel(rodata.srcsel ~= 0)).alpha];            % B1950
        decsrc = [rodata.catalog(rodata.srcsel(rodata.srcsel ~= 0)).delta];           % B1950
        epoch = true(length(rodata.srcsel), 1);
        % disp (['Sun NOT included, len of epoch ' num2str(length(epoch))]);
    end
    % srcpos0 = radectoITRF(rasrc(sel), decsrc(sel), epoch(sel), tobs);

    % srcpos0 holds the 3CR positions of ALL model sources, at the time of 
    % observation, and in ITRF coordinates. These are used as initial estimates 
    % for WSF based position estimation. 
    srcpos0 = radectoITRF(rasrc, decsrc, epoch, tobs_jd); 

    % Determine which sources are above the local horizon (visible to us).
    up = srcpos0 * rodata.normal > 0;
	sel = sigmas1 > 0.01 & up';

    % NOTE: We only estimate source positions for sources with an apparent flux
    % larger than 1% of the apparent flux of Cas A. 'sel' holds this subset of
    % sources with enough power, and above the horizon
    nsrc = sum(sel);
    if (calim.debug > 0)
       disp (['WSF: Working with ' num2str(nsrc) ...
             ' selected sources, tobs(MJD): ' num2str(tobs)]);
    end
    disp(['Srcs from srclist above horizon: ' num2str(rodata.srcsel(sel))]);


    % Convert coordinates of selected sources from ITRF to elevation/azimuth.
    thsrc0 = asin(srcpos0(sel, 3));
    phisrc0 = atan2(srcpos0(sel, 2), srcpos0(sel, 1));

	% Generate reference wsf estimation with full convergence
	calim.opt = optimset ();
	[thsrc_cat(:, maxiter), phisrc_cat(:, maxiter), thsrc_wsf(:, maxiter), ...
		phisrc_wsf(:, maxiter), fval, out] = ... 
       		 wsf_srcpos (simsky_acc, cal1, freq, sigman, ... 
						nsrc, rodata, phisrc0, ...
					    thsrc0, sel, calim.opt, calim.debug);
	fprintf (2, 'Ref: iters/calls/funcval, msg: %d / %d / %.2f, %s\n', ...
			out.iterations, out.funcCount, fval, out.message);
	
	% Iterate on the same timeslice of data multiple times 
	for ind=1:maxiter-1 % Last position holds full convergence solution
		if (ind == 1) iters (ind) = 1;
		else iters (ind) = 5*ind; end;
		
		calim.opt = optimset ('MaxIter', iters(ind));
		% calim.opt = optimset (); % Uncomment to calibrate to convergence
		[thsrc_cat(:, ind), phisrc_cat(:, ind), thsrc_wsf(:, ind), ...
		phisrc_wsf(:, ind), fval, out] = ... 
       		 wsf_srcpos (simsky_acc, cal1, freq, sigman, ... 
						nsrc, rodata, phisrc0, ...
					    thsrc0, sel, calim.opt, calim.debug);
		wsf_output(:,ind) = [fval, out.funcCount, out.iterations];
		fprintf (1, '%d ', ind);
	end;

	for ind = 1:maxiter-1
		% Display outputs
		fprintf (1, '\nTrial: %d, Func. value: %f, calls: %d, Iters: %d\n', ... 
				 ind, wsf_output(1,ind), wsf_output(2,ind), wsf_output(3,ind)); 
		for src = 1:size(thsrc_cat,1)
			fprintf (1, 'src: %d, cat:(%7.4f, %7.4f), WSF:(%7.4f, %7.4f)\n', ...
		src, thsrc_cat(src, ind), phisrc_cat(src, ind), thsrc_wsf(src, ind), ...
			phisrc_wsf(src, ind)); 
		end;
		leg {ind} = sprintf ('iter:%02d', iters(ind));
	
		% plot residuals as a function of position.
		set (gca, 'Fontsize', 14);
		subplot (1,3,1)
%		plot (((thsrc_wsf (:,maxiter) - thsrc_wsf(:,ind))*180/pi), ...
%			  ((phisrc_wsf(:,maxiter) - phisrc_wsf(:,ind))* 180/pi), ...
%			  char (col(mod (ind, length(col)) + 1)));
		plot (hypot ((thsrc_wsf (:,maxiter) - thsrc_wsf(:,ind))*180/pi , ...
			   (phisrc_wsf(:,maxiter) - phisrc_wsf(:,ind))* 180/pi), ...
 	 			  char (col(mod (ind, length(col)) + 1)));
		xlabel ( 'Th fullconvwsf - Th wsf (deg)');
		ylabel ( 'Phi fullconvwsf - Phi wsf (deg)');
		hold on;
		legend (leg);

		subplot (1,3,2)
		plot (((thsrc_cat(:,ind) - thsrc_wsf(:,ind))*180/pi), ((phisrc_cat(:,ind) - phisrc_wsf(:,ind)) * 180/pi), char (col(mod (ind, length(col)) + 1)));
		xlabel ( 'Th cat - Th wsf (deg)');
		ylabel ( 'Phi cat - Phi wsf (deg)');
		hold on;

%		subplot (1,2,1);
%		plot (((thsrc_cat(:,ind) - thsrc_wsf(:,ind)) * 180/pi), char (col(mod (ind, length(col)) + 1)));
%		ylabel ('Theta (cat) - Theta wsf (deg)');
%		hold on;
%
%		subplot (1,2,2);
%		plot (((phisrc_cat(:,ind) - phisrc_wsf(:,ind)) * 180/pi), char (col(mod (ind, length(col)) + 1)));
%		ylabel ( 'Phi cat - Phi wsf (deg)');
	end;

	for ind = 1:size (thsrc_cat, 1)
		% Plot residuals as a function of iteration number.
		subplot (1,3,3)
		plot (((thsrc_wsf(ind, :) - thsrc_wsf(ind, maxiter))*180/pi), char (col(mod (ind, length(col)) + 1)));
% ((phisrc_cat(ind, :) - phisrc_wsf(ind, :)) * 180/pi), char (col(mod (ind, length(col)) + 1)));
		xlabel ( 'Iteration number');
		ylabel ( 'Phi fullconvwsf - Phi wsf (deg)');
		hold on;
	end;

	sol.tobs = tobs;
	sol.iters = iters;
	sol.freq = freq;
	sol.thsrc_wsf  = thsrc_wsf;
	sol.phisrc_wsf = phisrc_wsf;
	sol.thsrc_cat  = thsrc_cat;
	sol.phisrc_cat = phisrc_cat;
	sol.wsf_output =  wsf_output;
