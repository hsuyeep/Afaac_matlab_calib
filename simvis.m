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
% pep/28Apr15

function [out] = simvis (parm)
	if (isempty (parm))
		array_rad = 2.5; % meters
		Nelem2rad = 5;
		array_spacing = array_rad/Nelem2rad;
		deb = 1;
		l0 = 0.25, m0 = 0.25; % Locate sources of unit amplitude at the locations (in l,m units) as specified above.
		deb = 1;
		arrayconfig= 'rect';
	else
		array_rad = parm.array_rad;
		Nelem2rad = parm.Nelem2rad;
		array_spacing = array_rad/Nelem2rad;
		deb = parm.deb;
		l0 = parm.l0;
		m0 = parm.m0;
		arrayconfig = parm.arrayconfig;
	end;
	arraysampling = [-array_rad:array_spacing:array_rad]; % Array extent = -5:5m, elements of array placed every 0.5m (=lambda/2 for 1m lambda).

	fprintf (1, '<-- %s array config chosen.\n', arrayconfig);
	switch lower(arrayconfig)
		case 'rect'
			% Create a rectangular grid of antenna, equi-spaced
			[xpos, ypos] = meshgrid (arraysampling);

		case 'disk'
			% Create a disk of antennas, with radius being the extent of array sampling.
            [xpos, ypos] = meshgrid (arraysampling);
			discsel = (sqrt(xpos.^2 + ypos.^2) < array_rad);
			xpos = xpos (discsel);
			ypos = ypos(discsel);

		case 'ring'
			% Create a ring of antenna, with inner rad = inner_rad
            [xpos, ypos] = meshgrid (arraysampling);
			inner_rad = 3; % meters
			innerdiscsel = (sqrt(xpos.^2 + ypos.^2) < inner_rad);
			ringsel = discsel - innerdiscsel;
			xpos = xpos(ringsel==1);
			ypos = ypos(ringsel==1);

		case 'lba_outer'
			load ('poslocal_outer.mat', 'poslocal');
			xpos = poslocal(:,1) - poslocal(1,1);
			ypos = poslocal(:,2) - poslocal(1,2);

		case 'lba_inner'
			load ('poslocal_inner.mat', 'poslocal');
			xpos = poslocal(:,1) - poslocal(1,1);
			ypos = poslocal(:,2) - poslocal(1,2);

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


	% we = exp (2*pi*i*(xpos(1,:)'*l0 + ypos(:,1)*m0)); % Generate phasor due to location of each element
    we = exp (2*pi*i*(xpos(:)*l0 + ypos(:)*m0)); % Generate phasor due to location of each element
	V = we * we'; % Generate the visibilities for the system at hand.

	out.img_l = [-1:0.01:1];
	out.img_m = out.img_l;

	% Create an image using acm2skymap:
	% out.map = acm2skyimage (V, xpos(1,:)', ypos(:,1), 299792458, out.img_l, out.img_m);
    out.map = acm2skyimage (V, xpos(:), ypos(:), 299792458, out.img_l, out.img_m);
    out.V = V;
	if (deb > 0)
        figure();
        mask = meshgrid(out.img_l).^2 + meshgrid(out.img_m).'.^2 < 1;
		imagesc (out.img_l, out.img_m, real(out.map).*mask); colorbar; axis tight;
		xlabel ('l'); ylabel ('m');
		title ('Simulated map');
        
        figure();
        subplot (121);
        imagesc (abs(out.V)); colorbar; axis tight;
		xlabel ('Ant'); ylabel ('Ant');
		title ('ACM amplitude');
        subplot (122);
        imagesc (angle(out.V)); colorbar; axis tight;
		xlabel ('Ant'); ylabel ('Ant');
		title ('ACM phase');
	end;
	out.V = V;
