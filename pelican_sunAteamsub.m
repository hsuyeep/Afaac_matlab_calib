% Program to implement Sun and A-team subtraction, after carrying out a 
% StefCal stage to completion. The function operates on a single timeslice of a 
% single channel.
% pep/14Apr2012
% Added documentation.
% pep/26Sep12
% Changed return parameters to structures instead of scattered variables.
% pep/28Jan13.

% function [thsrc_cat, phisrc_cat, thsrc_wsf, phisrc_wsf, suncomps, calvis, ...
%           gainsol, sigmas, sigman, good] = pelican_sunAteamsub (acc, t_obs, ...
%           freq, uvflag, flagant, debug, ptSun)
function [currsol] = pelican_sunAteamsub (acc, t_obs, ...
          freq, uvflag, flagant, debug, ptSun, max_calext_iter, ...
		  max_gainsolv_iter, posfilename)

% Arguments:
%    acc    : Filled, square ACM of size NantxNant complex numbers, where 
%              Nant = total unflagged antennas.
%    t_obs  : Time in MJD secs, as extracted from a MeasurementSet.
%    freq   : The frequency of the observation in Hz.
%    uvflag : Mask of size NantXNant for flagged visibilities. 1 => ignore.
%    flagant: Numerical list of antennas which are flagged (1-reference), leading
%             to a reshaping of the ACM.
%    debug  : Level of debug messages: 
%             0 ==> Only very important messages
%             1 ==> More detail
%             2 ==> Even more detail
%    ptSun  : 0 ==> Carry out a Sparse reconstruction of the Solar model.
%             1 ==> Model the Sun as a point source.
% NOT USED.
%   calim   : Structure of parameters to the calibration. If supplied as an 
%			  empty variable, calim is internally generated.

% Return values:
%    thsrc_cat : Elevation (radians) of catalog positions of model sky sources, 
%                at t_obs. 
%    phisrc_cat: Azimuth   (radians) of catalog positions of model sky sources, 
%                at t_obs.
%    thsrc_wsf : Elevation (radians) of a model sky source, as estimated using 
%    			WSF on the data.
%    azisrc_wsf: Azimuth   (radians) of a model sky source, as estimated using 
%    			WSF on the data.
%    calvis    : Calibrated visibilities, ready for imaging
%    gainsol   : Estimated per-element complex gains
%    sigmas    : Estimated flux of model sources
%    sigman    : Estimated system temprature per baseline
%    good      : Boolean encoding goodness of calibration solution

    % Load the following matrices once only, not on every call.
    persistent first_call;
    persistent rodata;
    persistent calim;
	persistent uloc_flag;
	persistent vloc_flag;
	persistent modsky calsky subsky duv Nuv uvpad;

	% NOTE: Returning only the unflagged visibilities.
	% calvis = zeros (size (acc)); 
    if (isempty (first_call))
        fprintf (2, '   --- Non Tracking Stef. calibration to convergence ---   \n\n');
    	% ---- Initialize Read-only data ---- 
    	rodata.C       = 299792458;         % speed of light, m/s
        rodata.lon     = 6.869837540;       % longitude of CS002 in degrees
        rodata.lat     = 52.915122495;      % latitude of CS002 in degrees
        rodata.Nelem   = 288;               % Max. number of elements
    	% 3 x 1 normal vector to the station field
        rodata.normal  = [0.598753, 0.072099, 0.797682].'; % Normal to CS002

        disp ('Loading 3CR catalog and local antenna positions.');
        % load ('poslocal.mat', 'posITRF', 'poslocal'); 
		load (posfilename, 'posITRF', 'poslocal');
    	rodata.posITRF = posITRF;
    	rodata.poslocal = poslocal;
    	load srclist3CR;
    	rodata.catalog = srclist3CR;  		% Sky catalog to use.

        % A-team from 3CR catalog for defining the model sky
    	% Nsrc x 1 vector with source indices. -->NOTE<--: use 0 for the Sun
        rodata.srcsel =  [324, 283, 88, 179, 0]; 

    	% ---- Calibration and imaging parameters ---- 
	    calim.restriction    = 10;    % Avoid vis. below 'restriction' 
	                                  % wavelengths length.
	    calim.maxrestriction = 60;    % Avoid vis. above 'maxrestriction' 
	                                  % meters (NOT wavelengths!)
	    calim.debug = debug;      	  % Set debug level.
	    calim.rem_ants = rodata.Nelem;
	    calim.flagant = flagant;
		calim.opt = optimset ();     % Note: Using all default values!
 	
 	    % Calibration stopping conditions
   	    calim.diffstop = 1e-3;         % difference bet. calib. solutions
		if (isempty (max_calext_iter))
	    	calim.maxiter  = 10;            % Max. cal_ext iterations.
		else
			calim.maxiter = max_calext_iter;
		end;

		% Max. gainsolv. iterations.
		if (isempty (max_gainsolv_iter))
	    	calim.maxiter_gainsolv = 30;    
		else
	    	calim.maxiter_gainsolv = max_gainsolv_iter;
		end;
	   

    	clear srclist3CR;
    	clear posITRF, poslocal;			% Couldn't think of a better way

		if (debug > 3)
	    	uloc = meshgrid (rodata.poslocal(:,1)) - ... 
					meshgrid (rodata.poslocal (:,1)).';
	    	vloc = meshgrid (rodata.poslocal(:,2)) - ... 
					meshgrid (rodata.poslocal (:,2)).';
		    [uloc_flag, vloc_flag] = gen_flagged_uvloc (uloc, vloc, flagant); 
	    	duv = 2.5;
	    	Nuv = 500; %1000                % size of gridded visibility matrix
	    	uvpad = 512; %1024              % specifies if any padding needs to be added
			modsky = figure;
			calsky = figure;
			subsky = figure;
		end;
        first_call = 0;
    end

    % Only uncomment for debug with acc.mat
    % disp ('---> Loading acc.mat');
    % load ('acc.mat');

    % NOTE: Hardcoded!
    currsol.good = 1;                    % Optimistic approach!
	currsol.tobs = t_obs;
	currsol.freq = freq;
    
    % Basic sanity check.
    if (isempty (acc) || isempty (t_obs))
        disp ('ACM or time missing!');   
        return;
    end
    fprintf (1, 'pelican_sunAteamsub: t_obs: %.1f,  Freq: %03f\n', ...
			  t_obs, freq);

 	% ---- In case of flagged antennas, generate reshaped ACM ---- 
    if (length(flagant) ~= 0)   
      	disp (['Flagant: ' num2str(flagant)]);
	    calim.flagant = flagant;
    	antmask = zeros (size (acc));
    	posmask = zeros (size (rodata.posITRF));
    	calim.rem_ants = length(acc) - length(flagant);
    	for ind = 1:length(flagant)
    		antmask (flagant(ind), :) = 1; antmask (:,flagant(ind)) = 1;
    		posmask (flagant(ind), :) = 1;
    	end
    	acc = reshape (acc(antmask ~= 1), [calim.rem_ants, calim.rem_ants]);
    	rodata.posITRF_fl = reshape(rodata.posITRF(posmask ~=1), ... 
        							[calim.rem_ants, 3]);
    	calim.uvflag  = reshape (uvflag(antmask ~= 1), ... 
        				   [calim.rem_ants, calim.rem_ants]);
    	disp (['NOTE: ACM resized after flagging to ' num2str(size (acc))]);
		% fprintf (1, 'Rank/rcond of ACM: %.4f/%.4f\n', rank(acc), rcond(acc));	

		% ---- Generate flagged uv coordinates. ----
		u=meshgrid(rodata.posITRF_fl(:, 1))-meshgrid(rodata.posITRF_fl(:, 1)).';
	   	v=meshgrid(rodata.posITRF_fl(:, 2))-meshgrid(rodata.posITRF_fl(:, 2)).';
	    w=meshgrid(rodata.posITRF_fl(:, 3))-meshgrid(rodata.posITRF_fl(:, 3)).';
	    uvw = [u(:), v(:), w(:)];
	    calim.uvdist = sqrt(sum(uvw.^2, 2) - (uvw * rodata.normal).^2);
    else
      	rodata.posITRF_fl = rodata.posITRF;
		calim.uvflag = uvflag;
		% ---- Generate flagged uv coordinates. ----
		u=meshgrid(rodata.posITRF_fl(:, 1))-meshgrid(rodata.posITRF_fl(:, 1)).';
	   	v=meshgrid(rodata.posITRF_fl(:, 2))-meshgrid(rodata.posITRF_fl(:, 2)).';
	    w=meshgrid(rodata.posITRF_fl(:, 3))-meshgrid(rodata.posITRF_fl(:, 3)).';
	    uvw = [u(:), v(:), w(:)];
	    calim.uvdist = sqrt(sum(uvw.^2, 2) - (uvw * rodata.normal).^2);
    end


	t_obs_mjdsec = t_obs; % NOTE: Just so debug images can update with mjd sec.
    t_obs = t_obs/86400 + 2400000.5; % Convert from MJD secs. to JD day units

    
    % whitening of the array covariance matrix (required for DOA estimation)
    % tic
    acc = acc ./ sqrt(diag(acc) * diag(acc).');

    % Initial calibration
    [cal1, sigmas1, Sigman1] = statcal_stefcal (acc, t_obs, freq, ... 
        rodata, calim, calim.uvflag);

    % [cal1, sigmas1, Sigman1] = statcal_vlss(acc, t_obs, freq, posITRF, srcsel,
    %                          normal, 10, maxrestriction, eye(Nelem), catalog);

    if (calim.debug > 1)
    	disp('First round of antenna calibration completed, sigmas1:'); 
    	disp (sigmas1');
    end
    % toc

    
    % Actual source positions may differ from the known positions due to
    % ionospheric image distortions, so we use a high-resolution DOA estimation
    % technique (weighted subspace fitting, WSF) to estimate source positions
    % relative to the position of Cas A, which is assumed to be the brightest
    % source in the sky (CasA is also above the horizon all year, at LOFAR 
    % latitudes).

    % NOTE: We use source positions from 3C catalog as initial estimate
    % tic
    if (sum(rodata.srcsel == 0) ~= 0) % Include Solar model
        [raSun, decSun] = SunRaDec(t_obs);                          % J2000
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
    % srcpos0 = radectoITRF(rasrc(sel), decsrc(sel), epoch(sel), t_obs);

    % srcpos0 holds the 3CR positions of ALL model sources, at the time of 
    % observation, and in ITRF coordinates. These are used as initial estimates 
    % for WSF based position estimation. 
    srcpos0 = radectoITRF(rasrc, decsrc, epoch, t_obs); 

    % Determine which sources are above the local horizon (visible to us).
    up = srcpos0 * rodata.normal > 0;

    % NOTE: We only estimate source positions for sources with an apparent flux
    % larger than 1% of the apparent flux of Cas A. 'sel' holds this subset of
    % sources with enough power, and above the horizon
    sel = sigmas1 > 0.01 & up'; 
    nsrc = sum(sel);
    if (calim.debug > 0)
       disp (['WSF: Working with ' num2str(nsrc) ...
             ' selected sources, t_obs(MJD): ' num2str(t_obs)]);
    end
    disp(['pelican_sunAteasub.m: Srcs from srclist above horizon: ' num2str(rodata.srcsel(up))]);


    % Convert coordinates of selected sources from ITRF to elevation/azimuth.
    thsrc0 = asin(srcpos0(sel, 3));
    phisrc0 = atan2(srcpos0(sel, 2), srcpos0(sel, 1));
    
	% Estimate positions using WSF.
    [currsol.thsrc_cat, currsol.phisrc_cat, currsol.thsrc_wsf, ... 
		currsol.phisrc_wsf, fval, out] = ... 
    wsf_srcpos (acc, cal1, freq, Sigman1, nsrc, rodata, phisrc0, ... 
                thsrc0, sel, calim.opt, calim.debug); 

	% Convert WSF positions from azi/el to ITRF
    srcposhat = [cos(currsol.phisrc_wsf(sel)) .* cos(currsol.thsrc_wsf(sel)),...
                 sin(currsol.phisrc_wsf(sel)) .* cos(currsol.thsrc_wsf(sel)),...
                 sin(currsol.thsrc_wsf(sel))];
    
    % Use updated source positions to improve calibration
    % tic
    A = exp(-(2 * pi * 1i * freq / rodata.C)*(rodata.posITRF_fl * srcposhat.'));
    %[ghat, sigmas2, Sigman2] = cal_ext_stefcal(acc, A, sigmas1(sel), ...
    %							squeeze(abs(Sigman1)) > 0, calim);
    [sol, stefsol] = cal_ext_stefcal(acc, A, sigmas1(sel), ...
    							squeeze(abs(Sigman1)) > 0, calim);
	% Need to make a temporary storage in sol,to prevent cal_ext.. from 
	% overwriting previous solution contents..
	currsol.flagant = flagant;
	antmaskvec = zeros (1, rodata.Nelem);
	antmaskvec ([calim.flagant]) = 1;
	currsol.gainsol = zeros (1, rodata.Nelem);
	currsol.gainsol(antmaskvec == 0) = sol.gainsol;
	currsol.sigmas = zeros (length(rodata.srcsel), 1);
	currsol.sigmas (sel) = sol.sigmas;
	currsol.sigmas_statcal = zeros (length(rodata.srcsel), 1);
	currsol.sigmas_statcal = sigmas1;
	currsol.sigman = sol.sigman;
	currsol.calext_iters = sol.calext_iters;
	currsol.stefcal_iters = sol.stefcal_iters;
	currsol.stefcal_tol = sol.stefcal_tol;
	currsol.pinv_sol = sol.pinv_sol;
	currsol.stefsol = stefsol; % Store the gainsolv () solution parameters.

	% Added because using final gains for WSF estimation in tracking cal
	% leads to a ripple in the estimated gain phase differences.
	currsol.cal1 = zeros (1, rodata.Nelem);
	antmaskvec = zeros (1, rodata.Nelem);
	antmaskvec ([calim.flagant]) = 1;
	currsol.cal1 (antmaskvec == 0) = cal1;

    % New calibration solutions. (NOT NEEDED, cal2 = currsol.gainsol)
    % cal2 = (1 ./ ghat)';

    % Calibrate visibilities removing extended emission and system noise
    % acccal = (cal2' * cal2) .* (acc - squeeze(Sigman2));
    % acccal = (currsol.gainsol' * currsol.gainsol) .* (acc - squeeze(currsol.sigman));
    acccal = (sol.gainsol' * sol.gainsol) .* (acc - squeeze(currsol.sigman));
    % toc
    
    %% A-team subtraction
    % We start by removing the A-team sources. Since the Sun is a barely
    % resolved source given the AARTFAAC resolution of 0.8 degrees at 60 MHz,
    % it requires separate treatment.
    % The difference between the actual positions as estimated by WSF and the
    % positions in the catalog is sufficiently large to make subtraction based
    % on catalog positions fail, we use the source positions estimated above.
    % We reduce srcsel to sources that appear brighter than 1% of Cas A.
    %  Note that index value 0 is used to denote the Sun in srcsel.
    % tic

    % Actual source indices into srclist3CR, don't use directly!
    visibleAteamsun = rodata.srcsel(sel); 

    % construct model visibilities for A-team sources, choosing Point source 
    % model of the Sun, if required.
    if (ptSun == 1)
        disp ('pelican_sunAteamsub.m: Carrying out Point source subtraction of the Sun');

    	% sel holds the indices of all relevant sources
    	A = exp(-(2 * pi * 1i * freq / rodata.C) * ... 
				(rodata.posITRF_fl * srcposhat.')); 
    	% RAteam = A * diag(sigmas2) * A';
    	RAteam = A * diag(currsol.sigmas(sel)) * A';
    else
    	A = exp(-(2 * pi * 1i * freq / rodata.C) * ... 
                (rodata.posITRF_fl * srcposhat(visibleAteamsun ~= 0, :).'));
    	% RAteam = A * diag(sigmas2(visibleAteamsun ~= 0)) * A';
    	RAteam = A * diag(currsol.sigmas(visibleAteamsun ~= 0)) * A';
    end 

	% subtract A-team
	accsubAteam = acccal - RAteam;
	accsubSunS  = accsubAteam;

	if (debug > 3)
		disp ('--> Creating instantaneous sky images');
		% Image model sky: FOR DEBUG!
	    % flagant = [1:12, 51, 206]; 
	    % load ('poslocal.mat', 'posITRF', 'poslocal'); 
	    [radecmap, calmap, calvis, l, m] = ... 
			fft_imager_sjw_radec (RAteam (:), uloc_flag(:), vloc_flag(:), ... 
									duv, Nuv, uvpad, t_obs, freq, 0);
		mask = zeros (size (calmap));
		mask (meshgrid (l).^2 + meshgrid(m).'.^2 < 0.9) = 1;
		figure(modsky);
		imagesc (l, m, real(calmap .* mask));
		colorbar;
		title (sprintf ('Model sky before subtraction from calibrated visibilities: Time: %f', t_obs_mjdsec));
	
		figure(calsky);
	    [radecmap, acccalmap, calvis] = ... 
			fft_imager_sjw_radec (acccal(:), uloc_flag(:), vloc_flag(:), ... 
									duv, Nuv, uvpad, t_obs, freq, 0);
		imagesc (l, m, real (acccalmap.* mask));
		% imagesc (abs(accsubAteam));
		colorbar;

		title (sprintf ('Calibrated visibilities before model subtraction. Time: %f', t_obs_mjdsec));
	
	    [radecmap, calmap, calvis] = ... 
			fft_imager_sjw_radec (accsubAteam(:), uloc_flag(:), vloc_flag(:), ... 
									duv, Nuv, uvpad, t_obs, freq, 0);
		figure(subsky);
		imagesc (l,m,real(calmap.*mask));
		colorbar;
		title (sprintf ('Image after Sun + A team subtraction. Time: %f', ...
				t_obs_mjdsec));
		% [sel, sel_l, sel_m] = overplotcat (t_obs_mjdsec, rodata.catalog, ...
	    %								300, subsky, true);
		pause;
	end;

    if (calim.debug > 0)
    	fprintf (1, 'A-team subtracted. Final model source fluxes (Sigmas): ');
		disp (currsol.sigmas');
    end 
    % toc


    if (ptSun == 0)
      %% Sun subtraction
      % The Sun is modeled using a 3-sparse grid of pixels around the estimated
      % position of the Sun. 
      % tic
      phi0sun = atan2(srcposhat(visibleAteamsun == 0, 2), ... 
    				  srcposhat(visibleAteamsun == 0, 1));
      th0sun = asin(srcposhat(visibleAteamsun == 0, 3));
      % define grid of pixels
      delta = 0.01;
      res = 2 * delta / 8;
      [phisun, thsun] = meshgrid(phi0sun-delta:res:phi0sun+delta, ... 
    							  th0sun-delta:res:th0sun+delta);
      phith_dist = sqrt((phisun(:) - phi0sun).^2 + (thsun(:) - th0sun).^2);
      phisun = phisun(phith_dist < delta);
      thsun  = thsun(phith_dist < delta);
      sunpos = [cos(phisun) .* cos(thsun), ... 
    			sin(phisun) .* cos(thsun), ...
    		    sin(thsun)];
      
      % solar component estimation using sparse reconstruction
      A = exp(-(2 * pi * 1i * freq / rodata.C) * ... 
			 (rodata.posITRF_fl * sunpos.'));
      accsubAteam = double (accsubAteam);
      % fluxSun = 1.20 * sigmas2(visibleAteamsun == 0);
      fluxSun = 1.20 * currsol.sigmas(visibleAteamsun == 0);
      cvx_begin
          variable sigmaSun(size(A, 2));
          minimize norm(accsubAteam(:) - khatrirao(conj(A), A) * sigmaSun, 2);
          subject to
          norm(sigmaSun, 1) <= fluxSun;
      cvx_end
      
      % iterate, increasing resolution
      sigmaSun(sigmaSun < 1e-3) = 0;
      phi0sun = phisun(sigmaSun ~= 0);
      th0sun  = thsun(sigmaSun ~= 0);
      res2    = res / 3;
      [phigrid, thgrid] = meshgrid(-res2:res2:res2);
      phisun2 = kron(phi0sun, ones(9, 1)) +  ...
    			kron(ones(length(phi0sun), 1), phigrid(:));
      thsun2  = kron(th0sun, ones(9, 1)) + ...
    			kron(ones(length(th0sun), 1), thgrid(:));
      sunpos2 = [cos(phisun2) .* cos(thsun2), ...
    			 sin(phisun2) .* cos(thsun2), ...
    			 sin(thsun2)];

      A = exp(-(2 * pi * 1i * freq / rodata.C) * ... 
			(rodata.posITRF_fl * sunpos2.'));
      % fluxSun = 1.20 * sigmas2(visibleAteamsun == 0);
      fluxSun = 1.20 * currsol.sigmas(visibleAteamsun == 0);
      cvx_begin
          variable sigmaSun(size(A, 2));
          minimize norm(accsubAteam(:) - khatrirao(conj(A), A) * sigmaSun, 2);
          subject to
          norm(sigmaSun, 1) <= fluxSun;
      cvx_end
      
      % try subtraction
      sigmaSun(sigmaSun < 1e-3) = 0;
      RSunS = A * diag(sigmaSun) * A';
      accsubSunS = acccal - RAteam - RSunS;
      if (calim.debug > 0)
    	  disp('Sun subtracted using Sparse reconstruction.');
      end
      %toc

      % Storing both fluxes and positions as az/ele
      currsol.suncomps = zeros (size (sigmaSun, 1), 3); 
      currsol.suncomps (:, 1) = sigmaSun;
      currsol.suncomps (:, 2) = atan2 (sunpos2 (:,2), sunpos2 (:, 1));  % Azi.
      currsol.suncomps (:, 3) = acos (sunpos2 (:,3));
    else
      currsol.suncomps = zeros (1, 3);  % Placeholder for when a pt. src sun 
										% subtraction is adequate.
    end  % if (ptSun == 1) ... 

	% calvis (antmask == 0) = single(accsubSunS); % Convert to float
    currsol.calvis      = single(accsubSunS); % Convert to float
    % gainsol     = cal2;
%   currsol.sigmas		= zeros (size (sel));
%	currsol.sigmas(sel) = sol.sigmas;
%   currsol.sigman      = Sigman2;
