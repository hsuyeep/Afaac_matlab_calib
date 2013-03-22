% Script to generate test input for the cal_ext_stefcal () function, and to 
% compare results with input parameters.
% NOTE: The cal_ext_stefcal () and internal functional prototypes are changed
% from the previous version in order to pull out solution characteristics, so
% they may have to be reverted to the older versions.
% pep/13Feb13

function func_cal_ext ()

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
	rodata.srcsel =  [324, 283, 88, 179];  % NOTE: No Sun being used!!


	% ----- Initialize calibration and imaging parameters ----- 
    % calim.restriction    = 10; % Avoid vis. below 'restriction' 
    calim.restriction    = 0; % Avoid vis. below 'restriction' 
                               % wavelengths length.
    calim.maxrestriction = 60; % Avoid vis. above 'maxrestriction' 
                               % meters (NOT wavelengths!)
    calim.debug = 0;      	   % Set debug level.
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


	% ----- Input parameters for simulation ----- 
	tobs = 4848774150.503317; % MJD sec, taken from LBA_OUTER_BAND_SPREAD obs.
	t_int  = 1; % 1 sec. integration
	f_int = 90000; % 90 KHz integration
	tobs_jd = tobs/86400.+2400000.5; % Convert to JD.
	freq = 69912719.726562;  % Hz
	duv = rodata.C/(freq*2); % Always Nyquist sample.
	abovejy = 300; 			 % Flux cutoff for generating a simulated sky.
	flagant = calim.flagant; % Flagged antennas.

	% This defines the variance of the gaussian noise to be added to the 
	% simulated complex gains per antenna ON TOP OF THE STATION GAINS.
	% antenna. This gives each antenna a random gain around a mean station gain,
	% (due to stations being more or less calibrated), thus simulating the 
	% AARTFAAC situation.
	antgainsig  = 0.5*ones (1, rodata.Nelem);

	% 1. No antenna gains at all.
	% fprintf (1, 'No antenna gains applied.\n');
	% stationgains = [0*exp(0*i*0.3), 0*exp(0*i*2.4), 0*exp(0*i*-0.0),...
 	% 			    0*exp(0*i*-0.2), 0*exp(0*i*2.8), 0*exp(0*i*-0.8)];

	% 2. Diverse antenna gains
	 fprintf (1, 'Diverse antenna gains applied, diff. amp/ph per station.\n');
	 stationgains = [1.2*exp(1*i*1.3), 0.8*exp(1*i*2.4), 0.5*exp(1*i*-1.1),...
  				    1.5*exp(1*i*-1.2), 1.8*exp(1*i*2.8), 0.7*exp(1*i*-1.8)];

	% 3. Only amplitude in antenna gains, no phases.
	% fprintf (1, 'Only ampl. in antenna gains applied, no phases.\n');
	%stationgains = [1.2*exp(0*i*1.3), 0.8*exp(0*i*2.4), 0.5*exp(0*i*-1.1),...
  	%			    1.5*exp(0*i*-1.2), 1.8*exp(0*i*2.8), 0.7*exp(0*i*-1.8)];
	
	% 4. Constant 1 as amp. gain, no phases.
	% fprintf (1, 'Const. 1 as ampl. gains applied, no phases.\n');
	% stationgains = [1*exp(0*i*1.3), 1*exp(0*i*2.4), 1*exp(0*i*-1.1),...
	% 			    1*exp(0*i*-1.2), 1*exp(0*i*2.8), 1*exp(0*i*-1.8)];

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

	% In m^2/K at 60MHz, and in the direction of CasA, from 'In situ antenna 
	% performance evaluation of the LOFAR phased array telescope', 
	% Wijnholds et. al, 2010
	% NOTE: this sensitivity can change depending on the location of the
	% galactic center (in boresight => low sensitivity, below horizon => high
	% and the zenith angle (boresight = high, low elevation = low)
	LBA_Aovertsys = 2.5e-3; % m^2/K, per LBA dipole at 60MHz. 

	% In Jy, 2*k*Tsys/Aeff, after putting in integration of 1s and 90 KHz
	LBA_dipole_sefd = (2*1.38*1e-3/LBA_Aovertsys)/sqrt(2*9e4); 

	% Extract out catalog srcs for given tobs and abovejy, above the horizon.
	% NOTE: Using 3CR fluxes directly, no spectral index.
	sel = ([rodata.catalog(:).flux] > abovejy);
	model_sky_fluxes = [rodata.catalog(sel).flux];

	model_sky_snr = model_sky_fluxes .* (LBA_Aovertsys/1.38*sqrt(9e4/2));
	disp (['Model sky SNR:' num2str(model_sky_snr)]);

	% sigma, in Jy, on observed visibility due to T_sys
	model_sky_vis_std = sqrt(2)*1.38*1e6/LBA_Aovertsys; 

%
%	% Gains have no units, just multiplicative factors!
%	total_Ta = total_flux^2*(antgains*antgains'); 
%	SNR = total_Ta/(LBA_dipole_sefd*288);
%	
%	% Now generating noise such that SNR over normalized total flux is identical
%	sigman = (1/SNR) * eye (rodata.Nelem);

%{
	LBA_beam_solid_angle_sr = (120*120)/((180/pi)^2); % 120deg 6dB beam at 60MHz
	t_gal = 60; % *lambda^2.55, from Bregman's thesis, in K;
	gal_flux = sqrt(2) * 1.38*1e-3*t_gal*(rodata.C/freq)^0.55*LBA_beam_solid_angle_sr; % *288*sqrt(f_int*t_int);
%In Jy
	t_gal = 60*(rodata.C/freq)^2.55; % From bregman's thesis
	casA_flux = 10^(5.625-0.634*log10(freq/1e6 - 0.023*(log10(freq/1e6)^2)));
	snr_casA = (casA_flux*1e-3/(2*1.38))*(LBA_Aovertsys*288)*sqrt(f_int*t_int);
	snr_gal = (gal_flux*1e-3/(2*1.38))*(LBA_Aovertsys*288)*sqrt(f_int*t_int);
	snr_allsky = 1e-4; % Arbit number for now!
%}

	% Control system noise. NOTE: Noise covariance should be Hermitiean
	% symmetric!
	% 1. All elements have identical, but diagonal noise, with std. in Jy.
	sigman = model_sky_vis_std*eye (rodata.Nelem); 

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


	% ----- Generation of simulated sky -----
    [simsky_acc, simsky_A, mod_acc, mod_A] = ... 
		genmodelsky (tobs, freq, [], rodata, ...
			calim, sigman, antgains, abovejy, 4); % calim.debug);

%	% Extract out 3CR sources for given tobs and abovejy, and above the horizon.
%	sel = ([rodata.catalog(:).flux] > abovejy);
%	srcpos_simsky =  ...
%		radectoITRF([rodata.catalog(sel).alpha], [rodata.catalog(sel).delta],...
%							true(sum(sel), 1), tobs_jd);
%    simsky_up = srcpos_simsky * rodata.normal > 0;
%	fprintf (1, 'testcal_ext:Found %d sources above %d Jy in visible sky\n', ...
%			 sum(simsky_up), abovejy);
%	catfluxrat = [rodata.catalog(sel).flux];	% NOTE: No Sun in simulations!
%
%	% CasA is the brightest 3CR src.
%	catfluxrat = catfluxrat ./ max (catfluxrat);
%	
%	% Generate synthetic visibilities for the model sky, with per-antenna
%	% gains and system noise applied.
%	simsky_A = exp(-(2*pi*1*i*freq/rodata.C)*(posITRF*srcpos_simsky(simsky_up,:).'));
%	simsky_acc = diag(antgains) * simsky_A*diag(catfluxrat(simsky_up))* ...
%				simsky_A' * diag (antgains)' + sigman;
%    simsky_acc = reshape (simsky_acc (antmask ~= 1), [calim.rem_ants, calim.rem_ants]);
%    disp (['NOTE: ACM resized after flagging to ' num2str(size (simsky_acc))]);
	

	% Generate details and image of simulated data
	fprintf (1, 'Norm/rank/rcond number: %f %f %f\n', norm(simsky_acc), ...
			 rank (simsky_acc), rcond(simsky_acc));
	fprintf (1, 'Min/max, Re: (%f, %f), Im: (%f, %f)\n', ... 
			min (real(simsky_acc(:))), max(real(simsky_acc(:))), ...
			min (imag(simsky_acc(:))), max(imag(simsky_acc(:)))); 

	if (calim.debug > 0)
		figure;
		subplot (1, 2, 1);
		plot (real(simsky_acc(:)), imag(simsky_acc(:)), 'o');
		title ('Re Vs. Im');
		xlabel ('Real'); ylabel ('Imag');
		subplot (1, 2, 2);
		plot (calim.uvdist(:), abs(simsky_acc(:)), 'o');
		title ('Amp Vs. uvdist');
		xlabel ('Uvdist (m)'); ylabel ('Amp');
	
		figure;
		[radecmap, calmap, calvis] = ... 
			fft_imager_sjw_radec (simsky_acc (:), uloc(:), vloc(:), ... 
									duv, Nuv, uvpad, tobs, freq, 0);
		imagesc (20*log10(abs(calmap)));
		colorbar;
		title (sprintf ('Simulated sky with %d model sources(>%d Jy).\n', ...
				 sum(simsky_up), abovejy));
	end;
	

	% --- Now being computed within genmodelsky.m ---
	% ----- Generate model sky for calibration. ----- 
	catfluxrat = [rodata.catalog(rodata.srcsel).flux];
	catfluxrat = catfluxrat ./ max (catfluxrat);
    srcpos0 = radectoITRF([rodata.catalog(rodata.srcsel).alpha], ...
						  [rodata.catalog(rodata.srcsel).delta], ...
						  true (length(rodata.srcsel),1), tobs_jd); 
    mod_up = srcpos0 * rodata.normal > 0;
    mod_A = exp(-(2 * pi * 1i * freq / rodata.C)*(rodata.posITRF *...
				srcpos0(mod_up, :).'));
	mod_acc = mod_A * diag(catfluxrat (mod_up)) * mod_A';
	
	% Generate image of skymodel
	if (calim.debug > 0)
		figure;
		[radecmap, calmap, calvis] = ... 
			fft_imager_sjw_radec (mod_acc (:), uloc(:), vloc(:), ... 
									duv, Nuv, uvpad, tobs, freq, 0);
		imagesc (abs(calmap));
		colorbar;
		title (sprintf ('Calibration sky model with %d model sources.\n', ...
				 sum(mod_up)));
	end;

	% ----- Calibration of simulated data -----
	% Call statcal for the first round of calibration
    [cal1, sigmas1, Sigman1] = statcal_stefcal (simsky_acc, tobs_jd, freq, ... 
        rodata, calim, eye(rodata.Nelem));

    sel_mod = sigmas1 > 0.01 & mod_up'; 
	% call cal_ext_stefcal with simulated visibilities
    [sol, stefsol] = cal_ext_stefcal(simsky_acc, mod_A, sigmas1(sel_mod), ...
    							squeeze(abs(Sigman1)) > 0, calim);


	% ----- Comparison of estimated parameters with true ones ----- 
	% Plot comparison of estimated parameters versus true parameters.
	figure;
	subplot (4, 1, 1);
	% plot (1./abs(cal1), 'ro'); % statcal_stefcal () solution
	plot (1./abs(sol.gainsol), 'ro'); % cal_ext_stefcal () solution
	hold on;
	plot (abs(antgains), 'bo');% Actual gains.
	title ('Estimated gains (red) Vs. simulated gains (blue)');
	xlabel ('Antenna number');
	ylabel ('Gain');


	subplot (4, 1, 2);
	% plot (angle(cal1), 'ro'); % statcal_stefcal () solution
	plot (angle(sol.gainsol), 'ro'); % cal_ext_stefcal () solution
	% plot (angle(cal1), 'ro');
	hold on;
	plot (angle(antgains), 'bo');
	title ('Estimated phases (red) Vs. simulated phases (blue)');
	xlabel ('Antenna number');
	ylabel ('phase (rad)');

	% Plot parameters in comparison with true ones.
	subplot (4, 1, 3);
	estimated_vals = [abs(sol.gainsol) angle(sol.gainsol)];
	true_vals = [abs(antgains) angle(antgains)];
	% plot (estimated_vals./true_vals);
	plot (abs(sol.gainsol) - abs (antgains), '-bo'); 
	xlabel ('Antenna number');
	title ('Estimated - true gains');

	subplot (4, 1, 4);
	% plot (estimated_vals - true_vals);
	plot (angle (sol.gainsol) - angle (antgains), '-bo'); 
	xlabel ('Antenna number');
	ylabel ('Rad');
	title ('Estimated - true phases');

	pinv_sol = abs(pinv(true_vals) * estimated_vals - 1);
	fprintf (2, 'testcal_ext: range of pinv_sol = (%f, %f)\n',...
			 min (pinv_sol(:)), max(pinv_sol(:)));

	% Plot system noise
	figure;
	plot (abs(sol.sigman(:)), 'ro');
	hold on;
	plot (abs(sigman(:)), 'bo');
	title ('Estimated system noise (red) Vs. simulated system noise (blue)');
	xlabel ('Visibility number');
	ylabel ('System noise');
