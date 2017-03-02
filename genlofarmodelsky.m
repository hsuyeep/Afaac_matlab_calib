% Script to generate a model sky based on a LOFAR sky model, given a time and
% frequency specification.
% pep/20Feb17

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
% simsky_up: selection of model sources above the horizon at given time.

% function [simsky_acc, simsky_A,  simsky_up] = ... 
function genlofarmodelsky (t_obs, freq, nant, skymodfname, rodata, ...
					calim, sigman, antgains,  abovejy, debug)

    res = 10; % Hardcoded resolution (arcsec) for gaussian model building.
    extent = 4200; % Extent of model sources to be modelled.
    fprintf (1, 'test');
	if (isempty (t_obs))
		fprintf (1, 'genlofarmodelsky: No time provided, using current time.');	
        t_obs = datenum2mjdsec(now);
	end;

	if (isempty (freq))
		fprintf (1, 'genlofarmodelsky: No freq. provided, using 60MHz.');	
        freq = 60e6;
	end;

	if (isempty (nant))
		fprintf (1, 'genlofarmodelsky: No array information provided, assuming AARTFAAC-6.');	
        nant = 288;
	end;

	if (isempty (skymodfname))
		fprintf (1, 'genlofarmodelsky: No sky model has been provided, assuming A-Team_4.sky.\n');
        skymodfname = 'A-Team_4.sky';
	end;

	if (isempty (antgains))
		fprintf (1, 'genlofarmodelsky: No antenna gains provided, setting them to one.\n');
        antgains = diag (nant);
	end;

	if (isempty (sigman))
		fprintf (1, 'genlofarmodelsky: No system noise provided, setting them to one.\n');
        sigman = zeros (nant);
	end;

	% Lazy caller, lets fill it ourselves.
	if (isempty (rodata))
    	% ---- Initialize Read-only data ---- 
        disp ('genlofarmodelsky: Initializing rodata.');
    	rodata.C       = 299792458;         % speed of light, m/s
        rodata.lon     = 6.869837540;       % longitude of CS002 in degrees
        rodata.lat     = 52.915122495;      % latitude of CS002 in degrees
        rodata.Nelem   = nant;              % Max. number of elements
    	% 3 x 1 normal vector to the station field
        rodata.normal  = [0.598753, 0.072099, 0.797682].'; % Normal to CS002

        fprintf (1, 'genlofarmodelsky: Loading local antenna positions.\n');
        if (nant == 288)
            load ('poslocal.mat', 'posITRF', 'poslocal'); 
        elseif (nant == 576)
            load ('poslocal_afaac12_outer.mat', 'posITRF', 'poslocal'); 
        end;

    	rodata.posITRF = posITRF;
    	rodata.poslocal = poslocal;
        fprintf (1, 'genlofarmodelsky: Loading sky model from file %s.\n',...
                 skymodfname);
        [mod, rarng, decrng, modimg] = readlofarskymodel (skymodfname, freq,...
                                                res, extent, 0);
    	clear posITRF, poslocal;			% Couldn't think of a better way
	end;
	
    if (isempty (calim))
        calim.flagant = [];
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

    simsky_acc = zeros (nant);

    fprintf (1, 'genlofarmodelsky: Creating vis for visible patches\n');
    for patch = 1:length (rarng) % For every patch in the model.
        pos  = radectoITRF (rarng{patch}, decrng{patch}, true, tobs_jd);
        simsky_up = pos * rodata.normal > 0;
        if (sum (simsky_up) == 0)
            fprintf (2, ...
             'genlofarmodelsky: %s patch not visible for time %d.\n', ...
                mod(patch).name{1}, t_obs);
            continue;
        end;

	    simsky_A = exp(-(2*pi*1i*freq/rodata.C)*(rodata.posITRF_fl*pos.'));
	        % simsky_acc = simsky_acc +  diag(antgains) * simsky_A*diag(flux(simsky_up))* ...
		    % 		simsky_A' * diag (antgains)' + sigman;
	    simsky_acc = simsky_acc +  simsky_A*modimg{patch}*simsky_A';
    
        if (debug > 4)
            figure ;
            imagesc (rarng{patch}, decrng{patch}, modimg{patch}); colorbar;
            xlabel ('DEC (rad)');
            ylabel ('RA  (rad)');
            title (sprintf ('Model of %s', mod(patch).name{1}));
        end;

        [l,m] = radectolm (mod(patch).meanra, mod(patch).meandec, tobs_jd, rodata.lon, ...
                                rodata.lat, false);
       fprintf (1, 'genlofarmodelsky: %s at (ra/dec): [ %7.4f/%7.4f], (l/m): [%7.4f/%7.4f] accumulated.\n', ...
                mod(patch).name{1}, mod(patch).meanra, mod(patch).meandec, l, m);
    end;


	if (debug > 3)
    	load srclist3CR;
    	rodata.catalog = srclist3CR;  		% Sky catalog to use.
       	rodata.srcsel =  [324, 283, 88, 179]; 
        rasrc  = [rodata.catalog(rodata.srcsel(rodata.srcsel ~= 0)).alpha];      
        decsrc = [rodata.catalog(rodata.srcsel(rodata.srcsel ~= 0)).delta];      
        modflux= [rodata.catalog(rodata.srcsel(rodata.srcsel ~= 0)).flux];      
        epoch = true(length(rodata.srcsel), 1);
        srcpos0 = radectoITRF(rasrc, decsrc, epoch, tobs_jd); 
        up = srcpos0 * rodata.normal > 0.1;
        [l3c, m3c] = radectolm (rasrc, decsrc, tobs_jd, rodata.lon*180/pi, ...
                            rodata.lat*180/pi, true);
        fprintf (1, '<-- A-team 3C positions: \n');
        for ind=1:length(rasrc)
            fprintf (1, '    <-- %s: %.5f, %.4f, l,m: %.4f, %.4f\n', ...
            rodata.catalog(rodata.srcsel(ind)).name, rasrc(ind), decsrc(ind),...
            l3c(ind), m3c(ind)); 
        end;
    	A = exp(-(2 * pi * 1i * freq / rodata.C) * ... 
				(rodata.posITRF_fl * srcpos0.')); 
    	RAteam = A * diag(modflux) * A';
		
		fprintf ('<-- genlofarmodelsky: Creating model sky images');
	    uloc = meshgrid (rodata.poslocal(:,1)) - ... 
				meshgrid (rodata.poslocal (:,1)).';
	    vloc = meshgrid (rodata.poslocal(:,2)) - ... 
				meshgrid (rodata.poslocal (:,2)).';
	    wloc = meshgrid (rodata.poslocal(:,3)) - ... 
				meshgrid (rodata.poslocal (:,3)).';
		[uloc_flag, vloc_flag, wloc_flag] = gen_flagged_uvloc (uloc, vloc, ...
                                                wloc, calim.flagant); 
		gparm.type = 'pillbox';
		gparm.duv = 0.5; 
		gparm.Nuv = 2000;
		gparm.uvpad = 2048; 
		gparm.fft = 1;
        gparm.lim = 0;
        l = linspace (-1, 1, gparm.uvpad);
	    [radecmap, calmap, calvis] = ... 
			fft_imager_sjw_radec (simsky_acc(:), uloc_flag(:), vloc_flag(:), ... 
								  gparm, [], [], t_obs, freq, 0);
	    [radecmap_3c, calmap_3c, calvis_3c] = ... 
			fft_imager_sjw_radec (RAteam(:), uloc_flag(:), vloc_flag(:), ... 
								  gparm, [], [], t_obs, freq, 0);
		figure;
		imagesc (l, l, abs(calmap));
		colorbar;
		title (sprintf ('Recovered Simulated sky. %s, Freq: %.2f', ... 
			   datestr (mjdsec2datenum (t_obs_mjdsec)), freq));
	    axis equal
   		axis tight
   		set (gca, 'YDir', 'Normal'); % To match orientation with station images
   		set (gca, 'XDir', 'Reverse'); % To match orientation with station images

    	ylabel('South $\leftarrow$ m $\rightarrow$ North', 'interpreter', 'latex');
    	xlabel('East $\leftarrow$ l $\rightarrow$ West', 'interpreter', 'latex');
		figure;
		imagesc (l, l, abs(calmap_3c));
		colorbar;
		title (sprintf ('Recovered 3C Simulated sky. %s, Freq: %.2f', ... 
			   datestr (mjdsec2datenum (t_obs_mjdsec)), freq));
	    axis equal
   		axis tight
   		set (gca, 'YDir', 'Normal'); % To match orientation with station images
   		set (gca, 'XDir', 'Reverse'); % To match orientation with station images

    	ylabel('South $\leftarrow$ m $\rightarrow$ North', 'interpreter', 'latex');
    	xlabel('East $\leftarrow$ l $\rightarrow$ West', 'interpreter', 'latex');
    end;
end
