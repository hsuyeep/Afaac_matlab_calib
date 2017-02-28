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
    fprintf (1, 'test');
	if (isempty (t_obs))
		fprintf (1, 'genmodelsky: No time provided, using current time.');	
        t_obs = datenum2mjdsec(now);
	end;

	if (isempty (freq))
		fprintf (1, 'genmodelsky: No freq. provided, using 60MHz.');	
        freq = 60e6;
	end;

	if (isempty (nant))
		fprintf (1, 'genmodelsky: No array information provided, assuming AARTFAAC-6.');	
        nant = 288;
	end;

	if (isempty (skymodfname))
		fprintf (1, 'genmodelsky: No sky model has been provided, assuming A-Team_4.sky.\n');
        skymodfname = 'A-Team_4.sky';
	end;

	if (isempty (antgains))
		fprintf (1, 'genmodelsky: No antenna gains provided, setting them to one.\n');
        antgains = diag (nant);
	end;

	if (isempty (sigman))
		fprintf (1, 'genmodelsky: No system noise provided, setting them to one.\n');
        sigman = zeros (nant);
	end;

	% Lazy caller, lets fill it ourselves.
	if (isempty (rodata))
    	% ---- Initialize Read-only data ---- 
        disp ('genmodelsky: Initializing rodata.');
    	rodata.C       = 299792458;         % speed of light, m/s
        rodata.lon     = 6.869837540;       % longitude of CS002 in degrees
        rodata.lat     = 52.915122495;      % latitude of CS002 in degrees
        rodata.Nelem   = nant;              % Max. number of elements
    	% 3 x 1 normal vector to the station field
        rodata.normal  = [0.598753, 0.072099, 0.797682].'; % Normal to CS002

        fprintf (1, 'genmodelsky: Loading local antenna positions.\n');
        if (nant == 288)
            load ('poslocal.mat', 'posITRF', 'poslocal'); 
        elseif (nant == 576)
            load ('poslocal_afaac12_outer.mat', 'posITRF', 'poslocal'); 
        end;

    	rodata.posITRF = posITRF;
    	rodata.poslocal = poslocal;
        fprintf (1, 'genmodelsky: Loading sky model from file %s.\n',...
                 skymodfname);
        mod = readlofarskymodel (skymodfname, freq);
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

    

    % Create model visibilities for every patch in the skymodel. 
    for patch = 1:length(mod)
        fprintf (1, '<-- Working on modeling response to patch %s (%.0f components).\n', ...
mod(patch).name{1}, length(mod(patch).patch{1}));
        if (debug > 3)
            meanra  = mean (mod(patch).patch{17});
            meandec = mean (mod(patch).patch{18});
            extent  = 6400; % arcsec
            modimg = zeros (extent/res);
            nrows = size (modimg, 1);
            ncols = size (modimg, 2);
            [l,m] = radectolm (meanra, meandec, tobs_jd, rodata.lon*180/pi, ...
                                rodata.lat*180/pi, false);
            fprintf (1, '<-- Source ra/dec:%.4f, %.4f   l,m: %.4f, %.4f\n',meanra, meandec, l,m);
        end;

        tot_flux = 0;
        for jind = 1:length(mod(patch).patch{1}) % Get number of components.
            comp = mod(patch).patch;
            tmp = comp{2}(jind); % Can't use comp(2){jind} directly! :/

            % Generate the positions of all the points making up the source.
            if (strcmp (tmp{1},'POINT'))
                pos  = radectoITRF (comp{17}(jind), comp{18}(jind), true, tobs_jd);
                flux = comp{8}(jind);

            elseif (strcmp (tmp{1},'GAUSSIAN'))
                rasig  = (comp{14}(jind)/3600)*(pi/180);
                decsig = (comp{15}(jind)/3600)*(pi/180);
                if (rasig == 0 && decsig == 0)
                    rasig  = (1/3600)*(pi/180); % Set to 1 arcsec (arbit)
                    decsig = (1/3600)*(pi/180); % Set to 1 arcsec (arbit)
                elseif (rasig == 0)
                    rasig = decsig;
                elseif (decsig == 0)
                    decsig = rasig;
                end;
                
                [rapos, decpos, flux] = gengaussiansrc (comp{17}(jind), ...
                   comp{18}(jind), true, comp{8}(jind), rasig, decsig, ...
                   comp{16}(jind), (res/3600)*(pi/180));
                pos = radectoITRF ((rapos(:)), (decpos(:)), true, tobs_jd);
            end;

            tot_flux = tot_flux + comp{8}(jind);
            if (debug > 3)
                % Offsets are in pixel units.
                offra  = int32(((comp{17}(jind) - meanra)*180*3600/pi)/res);
                offdec = int32(((comp{18}(jind) - meandec)*180*3600/pi)/res);
                decstart = offdec+int32(nrows/2);
                rastart = offra+int32(ncols/2);
                modimg(decstart:decstart + size (flux,1)-1, ...
                       rastart :rastart  + size (flux,2)-1) = ...
                      modimg(decstart:decstart + size (flux,1)-1, ...
                             rastart :rastart  + size (flux,2)-1) + flux;
            end;
            
            simsky_up = pos * rodata.normal > 0;
            if (sum (simsky_up) == 0)
                fprintf (2, 'Not visible.');
                % fprintf (2, '<-- Patch %s, comp %d not above the horizon for time %d.\n', mod(patch).name, jind, tobs_jd);
                continue;
            end;

	        simsky_A = exp(-(2*pi*1i*freq/rodata.C)*(rodata.posITRF_fl*pos.'));
	        % simsky_acc = simsky_acc +  diag(antgains) * simsky_A*diag(flux(simsky_up))* ...
		    % 		simsky_A' * diag (antgains)' + sigman;
	        simsky_acc = simsky_acc +  simsky_A*diag(flux(simsky_up))* ...
				simsky_A';
            
            % fprintf (1, '%d grid(%d,%d)   ', jind, size(flux,1), size(flux,2));
            fprintf (1, '%d ', jind);
        end;
        if (debug > 4)
            figure ;
            imagesc (modimg); colorbar;
            title (sprintf ('Model of %s', mod(patch).name{1}));
        end;
        fprintf (1, '..Done. Total flux: %f.\n', tot_flux);
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
        [l3c, m3c] = radectolm (rasrc, decsrc, tobs_jd, rodata.lon*180/pi, rodata.lat*180/pi, true);
        fprintf (1, '<-- A-team 3C positions: \n');
        for ind=1:length(rasrc)
            fprintf (1, '    <-- %s: %.5f, %.4f, l,m: %.4f, %.4f\n', rodata.catalog(rodata.srcsel(ind)).name, rasrc(ind), decsrc(ind), l3c(ind), m3c(ind)); 
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
		[uloc_flag, vloc_flag, wloc_flag] = gen_flagged_uvloc (uloc, vloc, wloc, calim.flagant); 
		gparm.type = 'pillbox';
		gparm.duv = 0.5; 
		gparm.Nuv = 2000;
		gparm.uvpad = 2048; 
		gparm.fft = 1;
        gparm.lim = 0;
	    [radecmap, calmap, calvis] = ... 
			fft_imager_sjw_radec (simsky_acc(:), uloc_flag(:), vloc_flag(:), ... 
								  gparm, [], [], t_obs, freq, 0);
	    [radecmap_3c, calmap_3c, calvis_3c] = ... 
			fft_imager_sjw_radec (RAteam(:), uloc_flag(:), vloc_flag(:), ... 
								  gparm, [], [], t_obs, freq, 0);
		figure;
		imagesc (abs(calmap));
		colorbar;
		title (sprintf ('Recovered Simulated sky. Time: %.2f, Freq: %.2f', ... 
			   t_obs_mjdsec, freq));
		figure;
		imagesc (abs(calmap_3c));
		colorbar;
		title (sprintf ('Recovered 3C Simulated sky. Time: %.2f, Freq: %.2f', ... 
			   t_obs_mjdsec, freq));
    end;
end

% Convert a gaussian model source to positions and fluxes 
% Arguments:
%  ra_cen : RA of gaussian peak, rad.
% dec_cen : Dec of gaussian peak, rad.
% epoch   : True => J2000
% tobs_jd : Time of observation in JD
% pk_flux : Peak flux of component, Jy.
% sigr/d  : Gaussian width in RA and DEC, in radians.
% pa      : Position angle in deg, in radians.
% res     : Resolution at which the Gaussian is to be constructed, in radians.
%
% Returns:
%   ragrid,decgrid: The positions of all parts of the component in radians    
function [ragrid, decgrid, flux] = gengaussiansrc (ra_cen, dec_cen, epoch, pk_flux, sigr, sigd, pa, res)

    extent = 3; % Extent of gaussian in sigma units.
    % Create a grid of points with resolution res radians, to cover the gaussian
    rarng = ra_cen + [-extent*sigr:res:extent*sigr];
    decrng = dec_cen + [-extent*sigd:res:extent*sigd];
    [ragrid, decgrid] = meshgrid (rarng, decrng);
    flux = pk_flux * exp (-( (ragrid-ra_cen).^2/(2*sigr^2) + (decgrid-dec_cen).^2/(2*sigd^2) ));
    % return ragrid, decgrid, flux;
end
