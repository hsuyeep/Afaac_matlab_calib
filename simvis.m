% Script to simulate visibilities from a specific array configuration, and source locations.
% Arguments:
% 	parm.array_rad : Extent of the array to be simulated
%	parm.Nelem2rad : Number of antenna elements in the array
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
% pep/28Apr15

function [out] = simvis (parm)
    def_freq = 60e6;  % Hz.
    
	if (exist ('parm') == 0)
		array_rad = 2.5; % meters
		Nelem2rad = 5;
		array_spacing = array_rad/Nelem2rad;
		deb = 1;
		l0 = 0.25, m0 = 0.25; % Locate sources of unit amplitude at the locations (in l,m units) as specified above.
		deb = 1;
		arrayconfig= 'rect';
        parm.fft = 1;
        parm.freq = def_freq;
	else
		arrayconfig = parm.arrayconfig;
        if (isfield (parm, 'fft') == 0)
            parm.fft = 1;
        end;
        if (isfield (parm,'freq') == 0)
            parm.freq = def_freq;
        end;

        % Need these parameters only for theoretical arrays
        if (strcmp (lower(parm.arrayconfig), 'lba_outer') == 0 && strcmp(lower(parm.arrayconfig),'lba_inner') == 0 && strcmp(lower(parm.arrayconfig), 'lba_outer_12') == 0)
		    array_rad = parm.array_rad;
		    Nelem2rad = parm.Nelem2rad;
		    array_spacing = array_rad/Nelem2rad;
	        arraysampling = [-array_rad:array_spacing:array_rad]; % Array extent = -5:5m, elements of array placed every 0.5m (=lambda/2 for 1m lambda).
        end;
		deb = parm.deb;
		l0 = parm.l0;
		m0 = parm.m0;
	end;

	fprintf (1, '<-- %s array config chosen.\n', arrayconfig);
	switch lower(arrayconfig)
        case 'linear'
            [xpos, ypos] = meshgrid (arraysampling, zeros (1, length (arraysampling)));

		case 'rect'
            % Example parameter structure:
            % parm.arrayconfig='rect'; parm.l0=0.1; parm.m0=0.25;
            % parm.array_rad=100; parm.Nelem2rad=5; parm.deb=1;

			% Create a rectangular grid of antenna, equi-spaced
			[xpos, ypos] = meshgrid (arraysampling);

		case 'disk'
            % Example parameter structure:
            % parm.arrayconfig='disk'; parm.l0=0.1; parm.m0=0.25;
            % parm.array_rad=100; parm.Nelem2rad=5; parm.deb=1;

			% Create a disk of antennas, with radius being the extent of array sampling.
            [xpos, ypos] = meshgrid (arraysampling);
			discsel = (sqrt(xpos.^2 + ypos.^2) < array_rad);
			xpos = xpos (discsel);
			ypos = ypos(discsel);

		case 'ring'
            % Example parameter structure:
            % parm.arrayconfig='ring'; parm.l0=0.1; parm.m0=0.25;
            % parm.array_rad=100; parm.Nelem2rad=5; parm.deb=1; parm.inner_rad =
            % 70;

			% Create a ring of antenna, with inner rad = inner_rad
            [xpos, ypos] = meshgrid (arraysampling);
			discsel = (sqrt(xpos.^2 + ypos.^2) < array_rad);
            if (isfield (parm, 'inner_rad') == 0)
                parm.inner_rad = 30;
            end;
			innerdiscsel = (sqrt(xpos.^2 + ypos.^2) < inner_rad);
			ringsel = discsel - innerdiscsel;
			xpos = xpos(ringsel==1);
			ypos = ypos(ringsel==1);

		case 'lba_outer'
			load ('poslocal_outer.mat', 'poslocal');
			xpos = poslocal(:,1); % - poslocal(1,1);
			ypos = poslocal(:,2); % - poslocal(1,2);

		case 'lba_inner'
			load ('poslocal_inner.mat', 'poslocal');
			xpos = poslocal(:,1) - poslocal(1,1);
			ypos = poslocal(:,2) - poslocal(1,2);

		case 'lba_outer_12'
			load ('poslocal_outer_afaac12.mat', 'poslocal');
			xpos = poslocal(:,1); % - poslocal(1,1);
			ypos = poslocal(:,2); % - poslocal(1,2);

		otherwise
			error (sprintf ('### Config %s is unimplemented!\n', arrayconfig));
	end;


	Nelem = length (xpos(:));

	if (deb > 0)
		% Show the array layout. A different colored dot for each row of elements in the array.
        figure();
		plot (xpos(:), ypos(:), '.'); 
		xlabel ('xpos (m)'); ylabel ('ypos(m)');
		title (sprintf ('Array configuration %s.\n', arrayconfig));
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
		imagesc (out.img_l, out.img_m, 10*log10(real(out.map)).*mask); colorbar; axis tight;
		imagesc (out.img_l, out.img_m, (real(out.map)).*mask); colorbar; axis tight;
		xlabel ('l'); ylabel ('m');
		title ('Simulated map, Ampl (dB)');
        
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
