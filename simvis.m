% Script to simulate visibilities from a specific array configuration, and source locations.
% Arguments:
% 	parm.arrayrad : Extent of the array to be simulated, in meters.
%	parm.elemspace : Element spacing, in meters.
%   parm.zslope    : Max. z height across elements in the array, in meters.
%   parm.pa        : Position angle of the array
%   parm.stretch   : Factor by which to stretch one dimension wrt. other.
%   parm.flux      : Flux (in arbit units) of specified sources in the sky.
%	parm.l0,m0     : Location of point sources in the sky, in l,m coordinates
%   parm.deb       : Bool controlling the generation of debug plots
%	parm.arrayconfig: Distribution of antennas. 
%					Currently supported:
%					- 'rect' : Equispaced rectangular 2-D grid of antennas
%					- 'disk' : Equispaced anntennas within a disk
%					- 'ring' : Equispaced antennas within an annuli
%					- 'lba_outer': LOFAR's LBA_OUTER actual configuration
%					- 'lba_outer_12': 12-station LBA_OUTER actual configuration
%					- 'lba_inner': LOFAR's LBA_INNER actual configuration
%  parm.fft        : Parameter defining whether FFT imaging should be carried
%                    out.
%  parm.freq       : Frequency of simulation. Needed as array extent is in
%                    meters.
%  parm.flagant    : Indices in the antenna array to be flagged.
% pep/28Apr15

function [out, parm] = simvis (parm)
    def_freq = 60e6;  % Hz.
    
    % If no parameters are given, this is our default case.
	if (exist ('parm') == 0)
		parm.arrayrad = 25;  % meters
		parm.elemspace = 2.5; % Default lambda / 2;
        parm.zslope = 0;
		parm.deb       = 1;
         % Locate sources of unit amplitude at the locations (in l,m units) as specified above.
        parm.flux = [5, 6];
		parm.l0 = [0.2, 0.25];
        parm.m0 = [0.3, -0.45];
		parm.deb       = 1;
		parm.arrayconfig = 'rect';
        parm.fft       = 1;
        parm.freq      = def_freq;
        parm.pa        = 0;
        parm.stretch   = 1;
        parm.ateam     = 0; % Switch off A-team simulation.
        parm.snr       = 5; % WGN to add to the visibilities.
        parm.flagant   = []; % None flagged.
        parm.tobs      = now();
	else
        % Check if all required parameters are available, else put in defaults,
        % even if unsed.
        if (isfield (parm, 'arrayrad') == 0)
            parm.arrayrad = 25;
        end;

        if (isfield (parm, 'elemspace') == 0)
            parm.elemspace = 2.5;
        end;

        if (isfield (parm, 'zslope') == 0)
            parm.zslope = 0;
        end;

        if (isfield (parm, 'pa') == 0)
            parm.pa = 0;
        end;

        if (isfield (parm, 'stretch') == 0)
            parm.stretch = 1;
        end;

        if (isfield (parm, 'l0') == 0)
            parm.flux = [5, 6];
        end;

        if (isfield (parm, 'l0') == 0)
		    parm.l0 = [0.2, 0.25];
        end;

        if (isfield (parm, 'm0') == 0)
            parm.m0 = [0.3, -0.45];
        end;

        if (isfield (parm, 'deb') == 0)
            parm.deb = 1;
        end;

        if (isfield (parm, 'arrayconfig') == 0)
            parm.arrayconfig = 'rect';
        end;

        if (isfield (parm, 'fft') == 0)
            parm.fft = 1;
        end;

        if (isfield (parm,'freq') == 0)
            parm.freq = def_freq;
        end;
        
        if (isfield (parm,'ateam') == 0)
            parm.ateam = 0;
        end;
        
        if (isfield (parm,'snr') == 0)
            parm.snr = 10;
        end;

        if (isfield (parm,'flagant') == 0)
            parm.flagant = [];
        end;

        if (isfield (parm,'tobs') == 0)
            parm.tobs = now;
        end;

    end;

    % Need these parameters only for theoretical arrays
    if (strcmp (lower(parm.arrayconfig), 'lba_outer') == 0 && strcmp(lower(parm.arrayconfig),'lba_inner') == 0 && strcmp(lower(parm.arrayconfig), 'lba_outer_12') == 0)
	    arraysampling_x = [-parm.arrayrad:parm.elemspace:parm.arrayrad]; 
        arraysampling_y = [-parm.arrayrad*parm.stretch:parm.elemspace:parm.arrayrad*parm.stretch];
        % Slope is always along the X dimension only.
        arraysampling_z = linspace (0, parm.zslope, length (arraysampling_x));
        normal = [0,0,1].';
    else
        normal = [0.598753, 0.072099, 0.797682].'; % Normal to CS002
    end;

    %%  Put in the Ateam sources for the given time.
    out.tobs = parm.tobs;
    out.freq = parm.freq;
    fprintf (2, '<-- Simulation for timestamp %s, freq %f.\n', datestr(out.tobs), out.freq);

    if (parm.ateam == 1)
	    % load the 3CR catalog for positions and fluxes
	    load 'srclist3CR.mat';
	    ateam_ind =  [324, 283, 88, 179];
        ateam_name = {'Cas.A', 'Cyg.A', 'Vir.A', 'Tau.A', 'Sun'};
	    ateam_ra = [srclist3CR(ateam_ind).alpha]; % RA
	    ateam_dec= [srclist3CR(ateam_ind).delta]; % DEC
	    ateam_fl = [srclist3CR(ateam_ind).flux];  % Jy
	    epoch = true ([1, length(ateam_ind)]);
	

        % Add Solar flux to the ateam
        [ateam_ra(length(ateam_ind)+1), ateam_dec(length(ateam_ind)+1)] = SunRaDec (datenum2mjdsec(out.tobs)/86400. + 2400000.5);

        % Flux of the quiet Sun in Jy.
        % from  [http://www.astro.phys.ethz.ch/astro1/Users/benz/papers/LBReview_thermal.pdf]
        % Valid 30<x<350MHz, freq. in MHz
        ateam_fl(length(ateam_ind)+1) = (1.94) * (out.freq/1e6)^1.992; 
        epoch(length(ateam_ind)+1) = false;

        [l, m] = radectolm(ateam_ra, ateam_dec, JulianDay(out.tobs), 6.869837540, 52.915122495,  epoch);
        fprintf (2, '<-- Sun located at RA/DEC: %.4f, %.4f, [l,m]: %.4f,%.4f\n', ateam_ra(end), ateam_dec(end), l(end), m(end));
	    % Convert ra/dec from catalog to ITRF coordinates
	    srcpos = radectoITRF(ateam_ra, ateam_dec, epoch, JulianDay (out.tobs));
	    up = srcpos * normal > 0.1; % Sources visible at this time.
        parm.l0 = l(up);
        parm.m0 = m(up);
        parm.flux = ateam_fl(up);
        parm.srcname = ateam_name(up);
        
	    % A = exp(-(2 * pi * 1i * obs.freq / 299792458)*(posITRF* srcpos(up).'));
	    % RAteam = A * diag(ateam_fl(sel)) * A';
    end;
    
    
	deb = parm.deb;
    assert (length(parm.l0) == length (parm.m0));
    assert (length(parm.flux) == length (parm.l0));
    fprintf (1, '<-- Simulated source properties (l,m,Jy(frac.power),Name): \n');
    for ind = 1:length(parm.l0)
        fprintf (1, '   [%5.2f, %5.2f, %8.2f (%6.4f), %s]\n', parm.l0(ind), parm.m0(ind), parm.flux(ind), parm.flux(ind)/sum(parm.flux),parm.srcname{ind});
    end 
    fprintf (1,'\n');
	l0 = parm.l0;
	m0 = parm.m0;

    % Generate the rotation matrix (only in 2D) for the given position angle.
    rotmat = [cosd(parm.pa) -sind(parm.pa);
              sind(parm.pa)  cosd(parm.pa)];

	fprintf (1, '<-- %s array config chosen.\n', parm.arrayconfig);
	switch lower(parm.arrayconfig)
        case 'linear'
            xpos = arraysampling_x; 
            ypos = ones (1, length(xpos)); 
            zpos = zeros (size (xpos));

		case 'rect'
            % Example parameter structure:
            % parm.arrayconfig='rect'; parm.l0=0.1; parm.m0=0.25;
            % parm.arrayrad=100; parm.Nelem2rad=5; parm.deb=1;

			% Create a rectangular grid of antenna, equi-spaced
			[xpos, ypos] = meshgrid (arraysampling_x, arraysampling_y);
            zpos = meshgrid (arraysampling_z, zeros (size (arraysampling_x)));

		case 'disk'
            % Example parameter structure:
            % parm.arrayconfig='disk'; parm.l0=0.1; parm.m0=0.25;
            % parm.arrayrad=100; parm.Nelem2rad=5; parm.deb=1;

			% Create a disk of antennas, with radius being the extent of array sampling.
            [xpos, ypos] = meshgrid (arraysampling_x, arraysampling_y);
			discsel = (sqrt(xpos.^2 + ypos.^2./parm.stretch) < parm.arrayrad);
			xpos = xpos (discsel);
			ypos = ypos (discsel);
            zpos = zeros (1, length (xpos));

		case 'ring'
            % Example parameter structure:
            % parm.arrayconfig='ring'; parm.l0=0.1; parm.m0=0.25;
            % parm.arrayrad=100; parm.Nelem2rad=5; parm.deb=1; parm.inner_rad =
            % 70;

			% Create a ring of antenna, with inner rad = inner_rad
            [xpos, ypos] = meshgrid (arraysampling_x, arraysampling_y);
			discsel = (sqrt(xpos.^2 + ypos.^2./parm.stretch) < parm.arrayrad);
            if (isfield (parm, 'inner_rad') == 0)
                parm.inner_rad = parm.arrayrad/2;
            end;
			innerdiscsel = (sqrt(xpos.^2 + ypos.^2./parm.stretch) < parm.inner_rad);
			ringsel = discsel - innerdiscsel;
			xpos = xpos(ringsel==1);
			ypos = ypos(ringsel==1);
            zpos = zeros (1, length (xpos));

		case 'lba_outer'
			load ('poslocal_outer_cs07w1.mat', 'poslocal');
			% load ('poslocal_outer_cs02w0.5.mat', 'poslocal');
			% load ('poslocal_outer.mat', 'poslocal');
            arraysampling_x = poslocal(:,1);
            arraysampling_y = poslocal(:,2);
            remants = setdiff([1:size(poslocal,1)], parm.flagant);
			xpos = poslocal(:,1); % - poslocal(1,1);
			ypos = poslocal(:,2); % - poslocal(1,2);
            zpos = poslocal(:,3);

		case 'lba_inner'
			load ('poslocal_inner.mat', 'poslocal');
			xpos = poslocal(:,1) - poslocal(1,1);
			ypos = poslocal(:,2) - poslocal(1,2);
            zpos = poslocal(:,3);

		case 'lba_outer_12'
            load ('poslocal_afaac12_outer_0w.mat', 'poslocal');
			% load ('poslocal_afaac12_outer.mat', 'poslocal');
			xpos = poslocal(:,1); % - poslocal(1,1);
			ypos = poslocal(:,2); % - poslocal(1,2);
            zpos = poslocal(:,3);

		otherwise
			error (sprintf ('### Config %s is unimplemented!\n', parm.arrayconfig));
	end;

    

    % Rotate the configuration
    tmp = rotmat * [xpos(:) ypos(:)]';
    xpos = reshape (tmp (1,:), size (xpos));
    ypos = reshape (tmp (2,:), size (ypos));

	Nelem = length (xpos(:));

	if (deb > 0)
		% Show the array layout. A different colored dot for each row of elements in the array.
        figure();
		plot3 (xpos(:), ypos(:), zpos(:), '.'); 
		xlabel ('xpos (m)'); ylabel ('ypos(m)'); zlabel ('zpos(m)');
		title (sprintf ('Array configuration %s.\n', parm.arrayconfig));
	end;

    uloc = (meshgrid (xpos(:)) - meshgrid (xpos(:)).'); % U in m 
    vloc = (meshgrid (ypos(:)) - meshgrid (ypos(:)).'); % V in m
    wloc = (meshgrid (zpos(:)) - meshgrid (zpos(:)).'); % W in m
    uvdist = sqrt (uloc(:).^2 + vloc(:).^2 + wloc(:).^2);
    V = sum(repmat ((parm.flux), size(uloc(:)),1) .* exp (-(2*pi*1i*parm.freq/299792458)*(uloc(:)*l0 + vloc(:)*m0 + wloc(:)*(sqrt(1-l0.^2-m0.^2) - 1))), 2);
    % V = sum(repmat ((parm.flux), size(uloc(:)),1) .* exp (-(2*pi*1i*parm.freq/299792458)*(uloc(:)*l0 + vloc(:)*m0 )), 2);
    V = conj (reshape (V, [length(xpos(:)), length(ypos(:))]));

    % Add noise to the visibilities
    if (parm.snr ~= 0)
        V = V + (mean(real(V(:)))/parm.snr) * randn (size (V)) +i*(mean(imag(V(:)))/parm.snr) * randn (size (V));
    end;
    
    % Another approach of first generating phasors in position coordinates.
    % Generate phasor due to location of each element
    % we = repmat ((parm.flux), size (xpos(:)), 1) .* (exp (-(2*pi*1i*parm.freq/299792458) *(xpos(:)*l0 + ypos(:)*m0) + zpos(:)*sqrt(1-l0.^2-m0.^2))); 
	% V = we * we'; % Generate the visibilities for the system at hand.

    if (parm.fft == 0)

	    % Create an image using acm2skymap:
	    out.img_l = [-1:0.01:1];
    	out.img_m = out.img_l;
        mask = meshgrid(out.img_l).^2 + meshgrid(out.img_m).'.^2 < 1;
        out.map = acm2skyimage (V, xpos(:), ypos(:), zpos(:), out.freq, out.img_l, out.img_m);
        out.map = flipud(rot90 (out.map)); % To match the locations of the
                                           % simulated point sources.
        scalefac = length (V(:));
    else

	    % Create a map using FFT imaging.
		gparm.type = 'pillbox';
	    gparm.lim  = 0;
	    gparm.duv = 0.5;				% Default, reassigned from freq. of obs. to
										% image just the full Fov (-1<l<1)
	    gparm.Nuv = 1500;				% size of gridded visibility matrix
	    gparm.uvpad = 1536;				% specifies if any padding needs to be added
	    gparm.fft  = 1;

	    % uloc = meshgrid (poslocal(:,1)) - meshgrid (poslocal (:,1)).';
	    % vloc = meshgrid (poslocal(:,2)) - meshgrid (poslocal (:,2)).';
	   	[radecmap, out.map, out.gridvis, out.img_l, out.img_m] = ... 
			  fft_imager_sjw_radec (V(:), uloc(:), vloc(:), ... 
						gparm, [], [], out.tobs, out.freq, 0);
        scalefac = (length(uloc(:))/gparm.Nuv).^2; % Still get a factor 10 extra
                                                   % flux... pep/03Oct16
    end;
    out.map = out.map/scalefac;
	out.V = V;
    out.uvdist = uvdist;
    out.xpos = xpos;
    out.ypos = ypos;
    out.zpos = zpos;

	if (deb > 0)
        figure();
        mask = meshgrid(out.img_l).^2 + meshgrid(out.img_m).'.^2 < 1;
		% imagesc (out.img_l, out.img_m, 10*log10(real(out.map)).*mask); colorbar; axis tight;
		imagesc (out.img_l, out.img_m, (real(out.map)).*mask); colorbar; axis tight;
		xlabel ('l'); ylabel ('m');
		title (sprintf ('Simulated map for %s array, max bline %.2f lambda, Ampl (dB)', parm.arrayconfig, max(uvdist(:))/(299792458/parm.freq)));
        
%        figure();
%        subplot (121);
%        imagesc (abs(out.V)); colorbar; axis tight;
%		xlabel ('Ant'); ylabel ('Ant');
%		title ('ACM amplitude');
%        subplot (122);
%        imagesc (angle(out.V)); colorbar; axis tight;
%		xlabel ('Ant'); ylabel ('Ant');
%		title ('ACM phase');
	end;
