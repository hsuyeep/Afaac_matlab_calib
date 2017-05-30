% Wrapper function to provide the appropriate data to cal_ext
% NOTE: CHANGE IN FUNCTION CALL to return sigmas (instead of sigmahat), which 
% has a zero for sources not above the horizon.

% Arguments:
% 	acc            : Nelem x Nelem x Nch array covariance matrix
% 	t_obs          : Nch x 1 vector of observing times
% 	freq           : Nch x 1 vector of observing frequencies
% 	rodata.posITRF : Nelem x 3 matrix with ITRF positions of the antennas
% 	rodata.srcsel  : Nsrc x 1 vector with source indices, use 0 for the Sun
% 	rodata.normal  : 3 x 1 normal vector to the station field
% 	calim.restriction: relative baseline restriction in wavelength
% 	calim.maxrestriction: maximum absolute baseline restriction in meters
% 	uvflag         : NelemxNelem bool matrix to reject(=0) specific baselines.
%   intaper        : Nelem x Nelem matrix to taper out the shortest baselines
%					 during estimation of gains, sigmas and select them for
%					 sigma_n estimation.
%
% Return values:
% 	cal    : Nelem x 1 vector with complex valued calibration corrections
% 	sigmas : Nsrc x 1 vector with estimated apparent source powers
% 	Sigman : Nelem x Nelem estimated array covariance matrix
%
% SJW, 2009
% modified on 18 May 2011 by SJW to use ITRF coordinates

function [cal, sigmas, Sigman, flux] = statcal_stefcal(acc, t_obs, freq, ...
					 rodata, calim, uvflag, mod_ra, mod_de)
% parameter section
	Nsb = length(freq);
	srcsel = rodata.srcsel;
	nsrc = length(srcsel);
	
	% initialization
	cal    = zeros(Nsb, calim.rem_ants);
	sigmas = zeros(Nsb, nsrc);
	Sigman = zeros(Nsb, calim.rem_ants, calim.rem_ants);
	
	for idx = 1:length(freq)
	    % disp (['NOTE: Time expected in MJD!']);
	    % disp(['working on subband 'num2str(idx) ' of 'num2str(length(freq))]);
		% Generate Array response matrix to model sky.
	    if (sum(srcsel == 0) ~= 0)
			% J2000
	        [raSun, decSun] = SunRaDec(t_obs(idx));                 

			% B1950
	        rasrc = [rodata.catalog(srcsel(srcsel ~= 0)).alpha, raSun].';

			% B1950
	        decsrc= [rodata.catalog(srcsel(srcsel ~= 0)).delta, decSun].'; 

	        epoch = [true(length(srcsel(srcsel ~= 0)), 1); false];
	    else
	        rasrc = [rodata.catalog(srcsel).alpha];  % B1950
	        decsrc= [rodata.catalog(srcsel).delta]; % B1950
	        epoch = true(length(srcsel), 1);
	    end
	    % disp (['t_obs: ' num2str(t_obs)]);
        if (~isempty (mod_ra))
            if (~isempty(mod_de))
                for modind = 1:length(mod_de)
                    fprintf (2, '<-- Additional model sky component at RA/Dec: %.2f, %.2f\n', mod_ra(modind), mod_de(modind));
                end;
				
                % rasrc = [rasrc mod_ra'];
                % decsrc = [decsrc mod_de'];
				% If model components are given, assume they are the brightest
				% sources in the sky, and do not attempt to estimate 
				% parameters for the Ateam sources.
                rasrc = mod_ra';
                decsrc = mod_de';
            end;
        end;
        
	    srcpos = radectoITRF(rasrc, decsrc, epoch, t_obs);
	    
	    up = srcpos * rodata.normal > 0.1;
	    A = exp(-(2 * pi * 1i * freq(idx) / rodata.C) * ... 
				(rodata.posITRF_fl * srcpos(up, :).'));
	    Rhat = squeeze(acc(:, :, idx));
	
		% Choose vis. subset based on visibility cut-off criteria and flags
	    % mask = reshape(calim.uvdist, [calim.rem_ants, calim.rem_ants]) < ... 
	    % 		   min([calim.restriction * (rodata.C / freq(idx)), ... 
		%			calim.maxrestriction]);
	    % mask =(1- mask) | uvflag; % Masked visibilties are used for tsys est.
		
        % Changed by Peeyush, 06Oct16, to apply a taper via downweighting the visibilities.
        % mask = calim.intap_fl & (1-uvflag); % if uvflag == 1, ignore vis.
		mask = calim.outtap_fl .* calim.intap_fl .* (1-uvflag); % if uvflag == 1, ignore vis.
	
		% Model source flux estimation via least squares imaging.
		% Ignore visibilities in the inner taper.
	    flux = real(((abs(A' * A).^2) \ khatrirao(conj(A), A)') * (Rhat(:) ... 
					.* mask (:)));
	    % flux = flux / flux(1);
		
		% For debug: pep/15Jan13
		% abs(A' * A)
		% figure; plot (khatrirao(conj(A), A)');
		
	    % flux(flux < 0) = 0;
	    if (calim.debug > 0)
	      disp ('Statcal_stefcal: model src flux estimates from LSimaging: ');
	      disp (flux');
	    end
	
	    % implementation of WALS method
	    % estimate direction independent gains, apparent source powers and
	    % receiver noise powers assuming that the source locations are correct
	    % [ghat, sigmahat, Sigmanhat]=cal_ext_stefcal(Rhat, A, flux, mask, calim);
		if (calim.debug > 1)
			fprintf (1, 'statcal_stefcal: Carrying out cal_ext.\n');
		end;
	    [sol, stefsol] = cal_ext_stefcal(Rhat, A, flux, (1-calim.intap_fl), uvflag, calim, 1, up);
	    
%	    cal(idx, :) = conj(1./ghat);
%	    sigmas(idx, up) = sigmahat; % Returning sigmahat
%	    Sigman(idx, :, :) = Sigmanhat;    
	    cal(idx, :) = sol.gainsol;
	    sigmas(idx, up) = sol.sigmas; % Returning sigmahat
	    Sigman(idx, :, :) = sol.sigman;    
	end
