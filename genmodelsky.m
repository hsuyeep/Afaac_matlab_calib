% Script to generate a model sky from the specified catalog, given a time and
% freq. specification.
% Note: If estimates of model source fluxes and positions are not provided,
% values from the specified catalog are used.
% pep/05Feb13

% Arguments: 
%  t_obs : Time in MJD secs of the model sky epoch.
%  freq  : Frequency in Hz. of the model sky observation
% modpos : The positions (in ITRF) of the model sources.
% rodata : Structure containing the required ancilliary data.
% calim  : Structure containing flagging and calibration related parameters.
% sigman : Noise contribution to visibilities.
%antgains: Complex gains associated with each antenna.
% abovejy: Flux cutoff of the population of sources in the model sky.
% debug  : Debug level of the script.

% Returns:
% simsky_acc: simulated sky visibilities.
% simsky_A : Array response to sources in simulated sky.
% mod_acc: model visibilities.
% mod_A : Array response to valid (above horizon etc.) model sources.
% mod_up: selection of model sources above the horizon at given time.

function [simsky_acc, simsky_A, mod_acc, mod_A, mod_up] = ... 
		genmodelsky (t_obs, freq, modpos, rodata, ...
					calim, sigman, antgains,  abovejy, debug)
	if (isempty (t_obs))
		fprintf (1, 'genmodelsky: Inconsistent time %f!', t_obs);	
	end;

	% Lazy caller, lets fill it ourselves.
	if (isempty (rodata))
    	% ---- Initialize Read-only data ---- 
        disp ('genmodelsky: Initializing rodata.');
    	rodata.C       = 299792458;         % speed of light, m/s
        rodata.lon     = 6.869837540;       % longitude of CS002 in degrees
        rodata.lat     = 52.915122495;      % latitude of CS002 in degrees
        rodata.Nelem   = 288;               % Max. number of elements
    	% 3 x 1 normal vector to the station field
        rodata.normal  = [0.598753, 0.072099, 0.797682].'; % Normal to CS002

        disp ('genmodelsky: Loading 3CR catalog and local antenna positions.');
        load ('poslocal.mat', 'posITRF', 'poslocal'); 
    	rodata.posITRF = posITRF;
    	rodata.poslocal = poslocal;
    	load srclist3CR;
    	rodata.catalog = srclist3CR;  		% Sky catalog to use.

        % A-team from 3CR catalog for defining the model sky
    	% Nsrc x 1 vector with source indices. -->NOTE<--: use 0 for the Sun
        rodata.srcsel =  [324, 283, 88, 179, 0]; 
    	clear srclist3CR;
    	clear posITRF, poslocal;			% Couldn't think of a better way

	end;
	
	% Handle flagged antennas
    antmask = zeros (rodata.Nelem);
    posmask = zeros (size (rodata.posITRF));
    rem_ants = rodata.Nelem - length(calim.flagant);
    for ind = 1:length(calim.flagant)
    	antmask (calim.flagant(ind), :) = 1; antmask (:,calim.flagant(ind)) = 1;
    	posmask (calim.flagant(ind), :) = 1;
    end
    % acc = reshape (acc(antmask ~= 1), [calim.rem_ants, calim.rem_ants]);
    rodata.posITRF_fl = reshape(rodata.posITRF(posmask ~=1), ... 
    							[rem_ants, 3]);

	% Convert from MJD secs. to MJD day units
	t_obs_mjdsec = t_obs;
    tobs_jd = t_obs/86400. + 2400000.5; 

	% NOTE: IGNORING modpos, inclusion still TODO.
	if (~isempty(modpos))
		fprintf (2, 'genmodelsky: Ignoring specified model source postions.\n');
	end;

	% Extract out catalog sources for given tobs and abovejy, and above the 
	% horizon.
	sel = ([rodata.catalog(:).flux] > abovejy);
	srcpos_simsky =  ...
		radectoITRF([rodata.catalog(sel).alpha], [rodata.catalog(sel).delta],...
							true(sum(sel), 1), tobs_jd);
    simsky_up = srcpos_simsky * rodata.normal > 0;
	fprintf (1, 'genmodelsky:Found %d sources above %d Jy in visible sky\n', ...
			 sum(simsky_up), abovejy);
	catfluxrat = [rodata.catalog(sel).flux];	% NOTE: No Sun in simulations!

	% CasA is the brightest 3CR src.
	% NOTE: DONT NORMALIZE, EVERYTHING IS IN JY NOW.
	% catfluxrat = catfluxrat ./ max (catfluxrat);
	
	% Generate Array response matrix, based on the number of sources, their
	% positions, and the positions of the array.
	simsky_A = exp(-(2*pi*1*i*freq/rodata.C)*(rodata.posITRF*srcpos_simsky(simsky_up,:).'));

	% Generate synthetic visibilities for the model sky, with per-antenna.
	% gains and system noise applied. Generated visibilities are in Jy, and are
	% without any noise contribution.
	simsky_acc = diag(antgains) * simsky_A*diag(catfluxrat(simsky_up))* ...
				simsky_A' * diag (antgains)' + sigman;

	% Add gaussian noise to real and imaginary components of generated vis.	
	% simsky_acc = simsky_acc + sigman*complex (randn(size(simsky_acc)), ...
	%				randn(size(simsky_acc)));
    simsky_acc = reshape (simsky_acc (antmask ~= 1), [calim.rem_ants, calim.rem_ants]);
    disp (['NOTE: ACM resized after flagging to ' num2str(size (simsky_acc))]);

%  --- Prev. catalog based sky generation ---
%	% Generate current positions of the model sources from catalog.
%    if (sum(rodata.srcsel == 0) ~= 0) % Include Solar model
%        [raSun, decSun] = SunRaDec(t_obs);                          % J2000
%        rasrc = [rodata.catalog(rodata.srcsel(rodata.srcsel ~= 0)).alpha, ...
%    			 raSun].';   % B1950
%        decsrc = [rodata.catalog(rodata.srcsel(rodata.srcsel ~= 0)).delta, ...
%    			  decSun].'; % B1950
%        epoch = [true(length(rodata.srcsel(rodata.srcsel ~= 0)), 1); false];
%        % disp (['Sun included, len of epoch ' num2str(length(epoch))]);
%    else
%        rasrc = [rodata.catalog(rodata.srcsel(rodata.srcsel ~= 0)).alpha];            % B1950
%        decsrc = [rodata.catalog(rodata.srcsel(rodata.srcsel ~= 0)).delta];           % B1950
%        epoch = true(length(rodata.srcsel), 1);
%        % disp (['Sun NOT included, len of epoch ' num2str(length(epoch))]);
%    end
%
%    % srcpos0 holds the 3CR positions of ALL model sources, at the time of 
%    % observation, and in ITRF coordinates. 
%    srcpos = radectoITRF(rasrc, decsrc, epoch, t_obs); 
%
%	% Generate model source flux ratios from catalog, if not given by user.
%	if (isempty (sigmas))
%		catfluxrat = [rodata.catalog(rodata.srcsel(rodata.srcsel ~= 0)).flux];	
%
%		% CasA is the brightest 3CR src.
%		catfluxrat = catfluxrat ./ max (catfluxrat);
%		disp (['genmodelsky: Flux ratios of model sources: '...
%				 num2str(catfluxrat)]);
%	end;
%
%    % Determine which sources are above the local horizon (visible to us).
%    up = srcpos * rodata.normal > 0;

    nsrc = sum(simsky_up);
    if (debug > 0)
      	fprintf (1, 'genmodelsky: Found %d src, t_obs(MJD): %.2f', nsrc, t_obs);
		%disp(['genmodelsky: Srcs from srclist above horizon: ' ...
% num2str(rodata.srcsel(simsky_up))]);
    end

	% ---- Generate model sky using calibrators only --- 
	% Extract out calibrator sources for given tobs and above the horizon.
	% Ignore Sun;
	fprintf (2, 'genmodelsky: Ignoring Sun!\n');
	modsel = rodata.srcsel (rodata.srcsel ~= 0);
	catfluxrat = [rodata.catalog(modsel).flux];
	catfluxrat = catfluxrat ./ max (catfluxrat);
    srcpos0 = radectoITRF([rodata.catalog(modsel).alpha], ...
						  [rodata.catalog(modsel).delta], ...
						  true (length(modsel),1), tobs_jd); 
    mod_up = srcpos0 * rodata.normal > 0;
    mod_A = exp(-(2 * pi * 1i * freq / rodata.C)*(rodata.posITRF_fl *...
				srcpos0(mod_up, :).'));
	mod_acc = mod_A * diag(catfluxrat (mod_up)) * mod_A';

%	if (isempty(sigmas)) % No model fluxes provided
%	    modsky_A = exp(-(2 * pi * 1i * freq / rodata.C)*(rodata.posITRF_fl * ...
%srcpos_modsky (modsky_up, :).'));
%		modsky_acc = modsky_A * diag(catfluxrat (modsky_up)) * modsky_A';
%	else
%    	sel = sigmas > 0.01 & up';  % Only choose sources with apparent power > 1% CasA.
%	    A = exp(-(2 * pi * 1i * freq / rodata.C)*(rodata.posITRF_fl * srcpos_simsky (sel, :).'));
%		simsky_acc = A * diag(sigmas(sel)) * A';
%	end;

	% Generate sky images if required.
	if (debug > 3)
		
		disp ('genmodelsky: Creating model sky images');
	    uloc = meshgrid (rodata.poslocal(:,1)) - ... 
				meshgrid (rodata.poslocal (:,1)).';
	    vloc = meshgrid (rodata.poslocal(:,2)) - ... 
				meshgrid (rodata.poslocal (:,2)).';
		[uloc_flag, vloc_flag] = gen_flagged_uvloc (uloc, vloc, calim.flagant); 
		gparm.type = 'pillbox';
		gparm.duv = 2.5; 
		gparm.Nuv = 500;
		gparm.uvpad = 512; 
		gparm.fft = 1;
	    [radecmap, calmap, calvis] = ... 
			fft_imager_sjw_radec (simsky_acc (:), uloc_flag(:), vloc_flag(:), ... 
								  gparm, [], [], t_obs, freq, 0);
		figure;
		imagesc (abs(calmap));
		colorbar;
		title (sprintf ('Simulated sky. Time: %.2f, Freq: %.2f', ... 
			   t_obs_mjdsec, freq));
	end;
