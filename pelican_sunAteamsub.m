% Program to implement Sun and A-team subtraction, after carrying out a 
% StefCal stage to completion. The function operates on a single timeslice of a 
% single channel.
% pep/14Apr2012
% Added documentation.
% pep/26Sep12

function [thsrc_cat, phisrc_cat, thsrc_wsf, phisrc_wsf, suncomps, calvis, ...
          gainsol, sigmas, sigman, good] = pelican_sunAteamsub (acc, t_obs, ...
          freq, uvflag, flagant, debug, ptSun)

% function [calvis, gainsol, sigmas, sigman, good] = 
%           pelican_sunAteamsub (acc, t_obs, freq, uvflag, debug)

%                         Function Arguments:
% acc    : Filled, square ACM of size NantxNant complex numbers, where 
%          Nant = total unflagged antennas.
% t_obs  : Time in MJD secs, as extracted from a MeasurementSet.
% freq   : The frequency of the observation in Hz.
% uvflag : Mask of size NantXNant for flagged visibilities. 1 => ignore.
% flagant: Numerical list of antennas which are flagged (1-reference), leading 
%          to a reshaping of the ACM.
% debug  : Level of debug messages: 
%          0 ==> Only very important messages
%          1 ==> More detail
%          2 ==> Even more detail
% ptSun  : 0 ==> Carry out a Sparse reconstruction of the Solar model.
%          1 ==> Model the Sun as a point source.

%                         Return values:
% thsrc_cat : Elevation (radians) of catalog positions of model sky sources, 
%             at t_obs. 
% phisrc_cat: Azimuth   (radians) of catalog positions of model sky sources, 
%             at t_obs.
% thsrc_wsf : Elevation (radians) of a model sky source, as estimated using WSF
%             on the data.
% azisrc_wsf: Azimuth   (radians) of a model sky source, as estimated using WSF
%             on the data.
% calvis    : Calibrated visibilities, ready for imaging
% gainsol   : Estimated per-element complex gains
% sigmas    : Estimated flux of model sources
% sigman    : Estimated system temprature per baseline
% good      : Boolean encoding goodness of calibration solution

    % Load the following matrices once only, not on every call.
    persistent first_call;
    persistent srclist3CR; % The 3rd Cambridge catalog (Revised)
    persistent posITRF;    % AARTFAAC Antenna locations (ITRF coordinates)
    persistent poslocal;   % AARTFAAC Antenna locations (wrt. local horizon)

    if (isempty (first_call))
        disp ('   --- Non Tracking Stef. calibration to convergence ---   ');
        disp ('Loading 3CR catalog and local antenna positions.');
    	load srclist3CR
    	load ('poslocal.mat', 'posITRF', 'poslocal'); 
        first_call = 0;
    end

    % Only uncomment for debug with acc.mat
    % disp ('---> Loading acc.mat');
    % load ('acc.mat');

    C = 299792458;               % speed of light, m/s
    lon = 6.869837540;           % longitude of CS002 in degrees
    lat = 52.915122495;          % latitude of CS002 in degrees
    restriction = 4;             % Avoid visibilities below 'restriction' 
                                 % wavelengths length.
    maxrestriction = 60;         % Avoid visibilities above 'maxrestriction' 
                                 % meters (NOT wavelengths!)

    % NOTE: Hardcoded!
    Nelem = 288;                 
    nblines = Nelem * (Nelem + 1)/2; 
    normal = [0.598753, 0.072099, 0.797682].';

    % A-team from 3CR catalog for defining the model sky
    srcsel =  [324, 283, 88, 179, 0]; 
    good = 1;                    % Optimistic approach!
	
    % Basic sanity check.
    if (isempty (acc) | isempty (t_obs))
        disp ('ACM or time missing!');   
        return;
    end
    disp (['t_obs :' num2str(t_obs) ' Freq: ' num2str(freq)]);

	% In case of flagged antennas, generate reshaped ACM.
    if (length(flagant) ~= 0)   
		disp (['Flagant' num2str(flagant)]);
		antmask = zeros (size (acc));
		posmask = zeros (size (posITRF));
		rem_ants = length(acc) - length(flagant);
		for ind = 1:length(flagant)
			antmask (flagant(ind), :) = 1; antmask (:,flagant(ind)) = 1;
			posmask (flagant(ind), :) = 1;
		end
		acc = reshape (acc(antmask ~= 1), [rem_ants, rem_ants]);
		posITRF_fl = reshape (posITRF(posmask ~=1), [rem_ants, 3]);
		uvflag  = reshape (uvflag(antmask ~= 1), [rem_ants, rem_ants]);
		disp (['NOTE: ACM resized after flagging to ' num2str(size (acc))]);
   	else
   	 	posITRF_fl = posITRF;
   	end

   	t_obs  = t_obs/86400 + 2400000.5; % Convert from MJD secs. to MJD day units
    catalog = srclist3CR;             % TODO: catalog selection to be removed?
    
    % whitening of the array covariance matrix (required for DOA estimation)
    % tic
    acc = acc ./ sqrt(diag(acc) * diag(acc).');

    % Initial calibration
	% TODO: put in restriction/maxrestriction here?
    [cal1, sigmas1, Sigman1] = statcal_stefcal (acc, t_obs, freq, posITRF_fl,...
                                  srcsel, normal, 10, 60, uvflag, debug);

    % [cal1, sigmas1, Sigman1] = statcal_vlss(acc, t_obs, freq, posITRF, srcsel,
    %                          normal, 10, maxrestriction, eye(Nelem), catalog);

    if (debug > 1)
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
    if (sum(srcsel == 0) ~= 0) % Include Solar model
        [raSun, decSun] = SunRaDec(t_obs);                          % J2000
        rasrc = [srclist3CR(srcsel(srcsel ~= 0)).alpha, raSun].';   % B1950
        decsrc = [srclist3CR(srcsel(srcsel ~= 0)).delta, decSun].'; % B1950
        epoch = [true(length(srcsel(srcsel ~= 0)), 1); false];
        % disp (['Sun included, len of epoch ' num2str(length(epoch))]);
    else
        rasrc = [srclist3CR(srcsel(srcsel ~= 0)).alpha];            % B1950
        decsrc = [srclist3CR(srcsel(srcsel ~= 0)).delta];           % B1950
        epoch = true(length(srcsel), 1);
        % disp (['Sun NOT included, len of epoch ' num2str(length(epoch))]);
    end
    % srcpos0 = radectoITRF(rasrc(sel), decsrc(sel), epoch(sel), t_obs);

	% srcpos0 holds the 3CR positions of ALL model sources, at the time of 
	% observation, and in ITRF coordinates. These are used as initial estimates 
	% for WSF based position estimation. 
    srcpos0 = radectoITRF(rasrc, decsrc, epoch, t_obs); 

	% Determine which sources are above the local horizon (visible to us).
    up = srcpos0 * normal > 0;

    % NOTE: We only estimate source positions for sources with an apparent flux
    % larger than 1% of the apparent flux of Cas A. 'sel' holds this subset of
	% sources with enough power, and above the horizon
    sel = sigmas1 > 0.01 & up'; 
    nsrc = sum(sel);
    if (debug > 0)
       disp (['WSF: Working with ' num2str(nsrc) ...
             ' selected sources, t_obs(MJD): ' num2str(t_obs)]);
    end
    disp(['Srcs from srclist above horizon: ' num2str(srcsel(up))]);


	% Convert coordinates of selected sources from ITRF to elevation/azimuth.
    thsrc0 = asin(srcpos0(sel, 3));
    phisrc0 = atan2(srcpos0(sel, 2), srcpos0(sel, 1));
    

    [thsrc_cat, phisrc_cat, thsrc_wsf, phisrc_wsf] = ... 
          wsf_srcpos (acc, cal1, freq, Sigman1, nsrc, posITRF_fl, phisrc0, ... 
	                  thsrc0, sel, srcsel, debug); 

    srcposhat = [cos(phisrc_wsf(sel)) .* cos(thsrc_wsf(sel)), ... 
                 sin(phisrc_wsf(sel)) .* cos(thsrc_wsf(sel)), ... 
                 sin(thsrc_wsf(sel))];
    
    % Use updated source positions to improve calibration
    % tic
    A = exp(-(2 * pi * 1i * freq / C) * (posITRF_fl * srcposhat.'));
    [ghat, sigmas2, Sigman2] = cal_ext_stefcal(acc, A, sigmas1(sel), ...
								squeeze(abs(Sigman1)) > 0, debug);
    % New calibration solutions.
    cal2 = (1 ./ ghat)';

    % Calibrate visibilities removing extended emission and system noise
    acccal = (cal2' * cal2) .* (acc - squeeze(Sigman2));
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
    visibleAteamsun = srcsel(sel); 

    % construct model visibilities for A-team sources, choosing Point source 
	% model of the Sun, if required.
    if (ptSun == 1)
    	disp ('Carrying out Point source subtraction of the Sun');

		% sel holds the indices of all relevant sources
		A = exp(-(2 * pi * 1i * freq / C) * (posITRF_fl * srcposhat.')); 
		RAteam = A * diag(sigmas2) * A';
    else
		A = exp(-(2 * pi * 1i * freq / C) * ... 
                (posITRF_fl * srcposhat(visibleAteamsun ~= 0, :).'));
		RAteam = A * diag(sigmas2(visibleAteamsun ~= 0)) * A';
    end 

    % subtract A-team
    accsubAteam = acccal - RAteam;
    accsubSunS  = accsubAteam;
    if (debug > 0)
		disp('A-team subtracted');
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
      A = exp(-(2 * pi * 1i * freq / C) * (posITRF_fl * sunpos.'));
      A = double (A);
      accsubAteam = double (accsubAteam);
      fluxSun = 1.20 * sigmas2(visibleAteamsun == 0);
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

      A = exp(-(2 * pi * 1i * freq / C) * (posITRF_fl * sunpos2.'));
      fluxSun = 1.20 * sigmas2(visibleAteamsun == 0);
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
      if (debug > 0)
		  disp('Sun subtracted using Sparse reconstruction.');
      end
      %toc

	  % Storing both fluxes and positions as az/ele
      suncomps = zeros (size (sigmaSun, 1), 3); 
      suncomps (:, 1) = sigmaSun;
      suncomps (:, 2) = atan2 (sunpos2 (:,2), sunpos2 (:, 1));  % Azimuth
      suncomps (:, 3) = acos (sunpos2 (:,3));
    else
      suncomps = zeros (1, 3);  % Placeholder for when a pt. src sun subtraction
                                % is adequate.
    end  % if (ptSun == 1) ... 


    calvis      = single(accsubSunS); % Convert to float
    gainsol     = cal2;
    sigmas(sel) = sigmas2;
    sigman      = Sigman2;

    % Quality of calibration solutions, as determined by the phase solutions
%    flag_ant_stations = zeros (1,6);
%    % Which station does current flagged antenna belong
%    for ant = 1:length(flagant)
%		flag_station = int32(flagant(ant)/48) + 1; 
%		flag_ant_stations (flag_station) = flag_ant_stations(flag_station) + 1;
%    end 
%    disp (['Flag_ant_stations ' num2str(flag_ant_stations)]);

%    ant_ind = 1;
%    for station=1:6
%    	ph_var (station) = var (angle(gainsol (1+(station-1)*48:station*48)));
%%      ph_var (station) = var (angle(gainsol (ant_ind:ant_ind+47-flag_ant_stations(station))));
%%      ant_ind = ant_ind + 47 -flag_ant_stations(station) + 1;
%    end
%    if sum(ph_var > 1) ~= 0 
%    	good = 0;
%    	disp ('NOTE: Bad calib solu:');
%    	disp (ph_var');
%    else
%    	good = 1;
%    end 
