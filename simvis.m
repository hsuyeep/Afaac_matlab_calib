% Script to simulate visibilities from a specific array configuration, and source locations.
% Arguments:
% 	parm.arrayrad : Extent of the array to be simulated, in meters.
%	parm.elemspace : Element spacing, in meters.
%   parm.pa        : Position angle of the array
%   parm.stretch   : Factor by which to stretch one dimension wrt. other.
%	parm.l0,m0     : Location of point sources in the sky, in l,m coordinates
%   parm.deb       : Bool controlling the generation of debug plots
%	parm.arrayconfig: Distribution of antennas. 
%					Currently supported:
%					- 'rect' : Equispaced rectangular 2-D grid of antennas
%					- 'disk' : Equispaced anntennas within a disk
%					- 'ring' : Equispaced antennas within an annuli
%					- 'lba_outer': LOFAR's LBA_OUTER actual configuration
%					- 'lba_inner': LOFAR's LBA_INNER actual configuration
%  parm.fft        : Parameter defining whether FFT imaging should be carried
%                    out.
%  parm.freq       : Frequency of simulation. Needed as array extent is in
%                    meters.
% pep/28Apr15

function [out, parm] = simvis (parm)
    def_freq = 60e6;  % Hz.
    
    % If no parameters are given, this is our default case.
	if (exist ('parm') == 0)
		parm.arrayrad = 25;  % meters
		parm.elemspace = 2.5; % Default lambda / 2;
		parm.deb       = 1;
         % Locate sources of unit amplitude at the locations (in l,m units) as specified above.
		parm.l0 = 0.25, parm.m0 = 0.25;
		parm.deb       = 1;
		parm.arrayconfig = 'rect';
        parm.fft       = 1;
        parm.freq      = def_freq;
        parm.pa        = 0;
        parm.stretch   = 1;
	else
        % Check if all required parameters are available, else put in defaults,
        % even if unsed.
        if (isfield (parm, 'arrayrad') == 0)
            parm.arrayrad = 50;
        end;

        if (isfield (parm, 'elemspace') == 0)
            parm.elemspace = 2.5;
        end;

        if (isfield (parm, 'pa') == 0)
            parm.pa = 0;
        end;

        if (isfield (parm, 'stretch') == 0)
            parm.stretch = 1;
        end;

        if (isfield (parm, 'l0') == 0)
            parm.l0 = 0.25;
        end;

        if (isfield (parm, 'm0') == 0)
            parm.m0 = 0.25;
        end;

        if (isfield (parm, 'deb') == 0)
            parm.deb = 1;
        end;

        if (isfield (parm, 'array_config') == 0)
            parm.array_config = 'rect';
        end;

        if (isfield (parm, 'fft') == 0)
            parm.fft = 1;
        end;

        if (isfield (parm,'freq') == 0)
            parm.freq = def_freq;
        end;
        
    end;
    % Need these parameters only for theoretical arrays
    if (strcmp (lower(parm.arrayconfig), 'lba_outer') == 0 && strcmp(lower(parm.arrayconfig),'lba_inner') == 0 && strcmp(lower(parm.arrayconfig), 'lba_outer_12') == 0)
	    arraysampling_x = [-parm.arrayrad:parm.elemspace:parm.arrayrad]; 
        arraysampling_y = [-parm.arrayrad*parm.stretch:parm.elemspace:parm.arrayrad*parm.stretch];
    end;
	deb = parm.deb;
	l0 = parm.l0;
	m0 = parm.m0;

    % Generate the rotation matrix (only in 2D) for the given position angle.
    rotmat = [cosd(parm.pa) -sind(parm.pa);
              sind(parm.pa)  cosd(parm.pa)];

	fprintf (1, '<-- %s array config chosen.\n', parm.arrayconfig);
	switch lower(parm.arrayconfig)
        case 'linear'
            [xpos, ypos] = meshgrid (arraysampling_x, arraysampling_y);

		case 'rect'
            % Example parameter structure:
            % parm.arrayconfig='rect'; parm.l0=0.1; parm.m0=0.25;
            % parm.arrayrad=100; parm.Nelem2rad=5; parm.deb=1;

			% Create a rectangular grid of antenna, equi-spaced
			[xpos, ypos] = meshgrid (arraysampling_x, arraysampling_y);
            zpos = zeros (1, length (xpos));

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
			load ('poslocal_outer.mat', 'poslocal');
			xpos = poslocal(:,1); % - poslocal(1,1);
			ypos = poslocal(:,2); % - poslocal(1,2);
            zpos = poslocal(:,3);

		case 'lba_inner'
			load ('poslocal_inner.mat', 'poslocal');
			xpos = poslocal(:,1) - poslocal(1,1);
			ypos = poslocal(:,2) - poslocal(1,2);
            zpos = poslocal(:,3);

		case 'lba_outer_12'
			load ('poslocal_outer_afaac12.mat', 'poslocal');
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
		plot (xpos(:), ypos(:), '.'); 
		xlabel ('xpos (m)'); ylabel ('ypos(m)');
		title (sprintf ('Array configuration %s.\n', parm.arrayconfig));
	end;

    uloc = (meshgrid (xpos) - meshgrid (xpos).'); % U in m 
    vloc = (meshgrid (ypos) - meshgrid (ypos).'); % V in m
    wloc = (meshgrid (zpos) - meshgrid (zpos).'); % W in m
    % V = exp (-(2*pi*1i*parm.freq/299792458)*(uloc*l0 + vloc*m0));
    
    % Generate phasor due to location of each element
    we = exp (-(2*pi*1i*parm.freq/299792458) *(xpos(:)*l0 + ypos(:)*m0)); 
	V = we * we'; % Generate the visibilities for the system at hand.

	out.tobs = now();
	out.freq = parm.freq;
    if (parm.fft == 0)

	    % Create an image using acm2skymap:
	    out.img_l = [-1:0.0025:1];
    	out.img_m = out.img_l;
        out.map = acm2skyimage (V, ypos(:), xpos(:), 299792458, out.img_l, out.img_m);
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
    end;
	out.V = V;

	if (deb > 0)
        figure();
        mask = meshgrid(out.img_l).^2 + meshgrid(out.img_m).'.^2 < 1;
		%imagesc (out.img_l, out.img_m, 10*log10(real(out.map)).*mask); colorbar; axis tight;
		imagesc (out.img_l, out.img_m, (real(out.map)).*mask); colorbar; axis tight;
		xlabel ('l'); ylabel ('m');
		title (sprintf ('Simulated map for %s array, max bline %d lambda, Ampl (dB)', parm.arrayconfig, parm.arrayrad/(299792458/parm.freq)));
        
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
