% Script to subtract sources from visibilities, based on the local coordinate positions specified.
% pep/14Apr15
% Arguments:
%  ra   : RA of sources (radians)
%  dec  : DEC of sources (radians)
%  epoch: false for J2000 positions
%  acm  : ACM
%  t_obs: time in MJDsec.
%  freq : Freq. of observation in Hz.
% uvflag: Mask to flag out chosen visibilities.
% flagant: list of antennas to be considered as flagged.
% deb   : Flag controlling the generation of debugging info.

function [accsub] = srcsub (ra, dec, epoch, acm, t_obs, freq, uvflag, flagant, deb)
        normal  = [0.598753, 0.072099, 0.797682].'; % Normal to CS002
    	C       = 299792458;         % speed of light, m/s
		load ('poslocal_outer.mat', 'posITRF', 'poslocal');
    	antmask = zeros (size (acm));
    	posmask = zeros (size (posITRF));
    	rem_ants = length(acm) - length(flagant);
        goodant = setdiff ([1:288], flagant);

    	for ind = 1:length(flagant)
    		antmask (flagant(ind), :) = 1; antmask (:,flagant(ind)) = 1;
    		posmask (flagant(ind), :) = 1;
    	end
    	acc = reshape (acm(antmask ~= 1), [rem_ants, rem_ants]);
    	posITRF_fl = reshape(posITRF(posmask ~=1), [rem_ants, 3]);
		
		% Convert source positions to ITRF
        tobs_jd = t_obs./86400 + 2400000.5; % Convert to JD
	    srcpos = radectoITRF(ra, dec, epoch, tobs_jd);
	    
		% Check if sources are above the horizon
	    up = srcpos * normal > 0;
        if (sum(up) == 0)
            warning ('No source found above horizon!');
        end;

		% Generate array manifold
	    A = exp(-(2 * pi * 1i * freq / C) * ... 
				(posITRF_fl * srcpos(up, :).'));
		mask = (1-uvflag(goodant, goodant)); % if uvflag == 1, ignore vis.
		% Flux estimate
	    flux = real(((abs(A' * A).^2) \ khatrirao(conj(A), A)') * (acc(:) ... 
					.* mask (:)));

        flux (flux < 0) = 0;
		% Generate model visibilities with estimated flux.
    	modvis = A * diag(flux) * A';

		accsub= acc - modvis;

		if (deb > 0)
	    	uloc = meshgrid (poslocal(:,1)) - ... 
					meshgrid (poslocal (:,1)).';
	    	vloc = meshgrid (poslocal(:,2)) - ... 
					meshgrid (poslocal (:,2)).';
		    [uloc_flag, vloc_flag] = gen_flagged_uvloc (uloc, vloc, flagant); 

			gparm.type = 'pillbox';
			gparm.lim = 0;
			gparm.duv = 0.5; 
			gparm.Nuv = 500;
			gparm.uvpad = 512; 
			gparm.fft = 1;

			% Show the raw acm image
			[map, calmap, calvis, l, m] = ... 
                    fft_imager_sjw_radec (acc, uloc_flag(:), vloc_flag(:), ...
										gparm, [], [], t_obs, freq, 0);
	
            [modmap, calmap, calvis, l, m] = ... 
                    fft_imager_sjw_radec (modvis, uloc_flag(:), vloc_flag(:), ...
										gparm, [], [], t_obs, freq, 0);
	                                    
                                    
			[mapsub, calmap, calvis, l, m] = ... 
                    fft_imager_sjw_radec (accsub, uloc_flag(:), vloc_flag(:), ...
										gparm, [], [], t_obs, freq, 0);

			subplot (1,3,1);
			imagesc (l,m,abs(map)); colorbar;
    		set (gca, 'YDir', 'Normal'); % To match orientation with station images
		    set (gca, 'XDir', 'Reverse'); % To match orientation with station images
		    title (sprintf ('Input. Time: %s', datestr(mjdsec2datenum(t_obs))));
		    ylabel('South $\leftarrow$ m $\rightarrow$ North', 'interpreter', 'latex');
		    xlabel('East $\leftarrow$ l $\rightarrow$ West', 'interpreter', 'latex'); 

			subplot (1,3,2);
			imagesc (l,m,abs(modmap)); colorbar;
    		set (gca, 'YDir', 'Normal'); % To match orientation with station images
		    set (gca, 'XDir', 'Reverse'); % To match orientation with station images
		    title (sprintf ('Input. Time: %s', datestr(mjdsec2datenum(t_obs))));
		    ylabel('South $\leftarrow$ m $\rightarrow$ North', 'interpreter', 'latex');
		    xlabel('East $\leftarrow$ l $\rightarrow$ West', 'interpreter', 'latex'); 

            subplot (1,3,3);
			imagesc (l,m,abs(mapsub)); colorbar;
    		set (gca, 'YDir', 'Normal'); % To match orientation with station images
		    set (gca, 'XDir', 'Reverse'); % To match orientation with station images
		    title (sprintf ('Output. Time: %s',datestr(mjdsec2datenum(t_obs))));
		    ylabel('South $\leftarrow$ m $\rightarrow$ North', 'interpreter', 'latex');
		    xlabel('East $\leftarrow$ l $\rightarrow$ West', 'interpreter', 'latex'); 
		end;
