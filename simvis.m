% Script to simulate visibilities from a specific array configuration, and 
% source locations.
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
%  parm.tobs       : Simulated time of the observation (as a datenum)
%  parm.ateam      : Flag controlling simulation of a sky with the A-team+Sun
%  parm.lofarskymodel: Filename of the lofarsky model to use. If empty, a
%                    point source model will be used.
%  parm.snr        : Controls the SNR of the visibilities.
%  parm.gain       : Complex gains to be associated with the antennas.
%                    [nant x 1]
%  parm.sigman     : Complex coloured noise to be added to the visibilities.
%                    [nant x nant]
% pep/28Apr15

function [out, parm] = simvis (parm)
    def_freq = 60e6;  % Hz.
    
    % If no parameters are given, this is our default case.
	if (exist ('parm') == 0)
		parm.arrayrad = 25;  % meters
		parm.elemspace = 2.5; % Default lambda / 2;
        parm.zslope = 0;
        parm.pa        = 0;
        parm.stretch   = 1;
        parm.flux = [5, 6];

        % Locate sources of unit amplitude at the locations (in l,m units) as
        % specified above.
		parm.l0 = [0.2, 0.25];
        parm.m0 = [0.3, -0.45];
		parm.deb       = 1;
		parm.arrayconfig = 'rect';
        parm.fft       = 1;
        parm.freq      = def_freq;
        parm.flagant   = []; % None flagged.
        parm.tobs      = now(); % datenum
        parm.ateam     = 0;  % Switch off A-team simulation.
        parm.lofarskymodel = [];
        parm.snr       = 5;  % WGN to add to the visibilities.
        parm.gain      = [];
        parm.sigman    = [];
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
            parm.ateam = 1;
            parm.lofarskymodel = 'A-Team_4.sky'; % Use specified sky model.
        end;
        
        if (isfield (parm,'lofarskymodel') == 0)
            parm.lofarskymodel = [];
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

        if (isfield (parm,'gain') == 0)
            parm.gain = [];
        end;

        if (isfield (parm,'sigman') == 0)
            parm.sigman = [];
        end;
    end;
    out.tobs = parm.tobs;
    out.freq = parm.freq;
    normal = [0.598753, 0.072099, 0.797682].'; % Normal to CS002

    % ---------------------------------------------------------------%
    % ----------------------- Create model sky ----------------------%
    % ---------------------------------------------------------------%
    %%  Put in the Ateam sources for the given time.
    if (parm.ateam == 1)

        if (isempty (parm.lofarskymodel))
    	    % load the 3CR catalog for positions and fluxes
    	    load 'srclist3CR.mat';
    	    ateam_ind =  [324, 283, 179, 88];
            nsrc      = length (ateam_ind);
            ateam_name = {'Cas.A', 'Cyg.A', 'Vir.A', 'Tau.A', 'Sun'};
    	    ateam_ra = [srclist3CR(ateam_ind).alpha]; % RA
    	    ateam_dec= [srclist3CR(ateam_ind).delta]; % DEC
    	    ateam_fl = [srclist3CR(ateam_ind).flux];  % Jy
    	    epoch = true ([1, length(ateam_ind)]);
        else
            [mod, rarng, decrng, modimg] = readlofarskymodel ...
                            (parm.lofarskymodel, parm.freq, 10, 4200, 0);
            nsrc = length(mod);
            for ind = 1:length (mod)
                ateam_name {ind} = mod(ind).name{1};
                ateam_ra (ind)   = mod(ind).meanra;
                ateam_dec (ind)  = mod(ind).meandec;
                ateam_fl (ind)   = mod(ind).meanflux;
            end;
            ateam_name{nsrc+1} = 'Sun';
        end;

        tobs_jd = datenum2mjdsec (out.tobs)/86400. + 2400000.5;
	

        % Add Solar flux to the ateam
        [ateam_ra(nsrc+1), ateam_dec(nsrc+1)] = SunRaDec (tobs_jd);

        % Flux of the quiet Sun in Jy.
        % from  [http://www.astro.phys.ethz.ch/astro1/Users/benz/papers/\
        % LBReview_thermal.pdf]
        % Valid 30<x<350MHz, freq. in MHz
        ateam_fl(nsrc+1) = (1.94) * (out.freq/1e6)^1.992; 
        epoch(nsrc+1) = false;

        [l, m] = radectolm(ateam_ra, ateam_dec, tobs_jd, 6.869837540, ...
                            52.915122495,  epoch);
        fprintf (2, ...
            '<-- Sun located at RA/DEC: %.4f, %.4f, [l,m]: %.4f,%.4f\n',...
            ateam_ra(end), ateam_dec(end), l(end), m(end));

	    % Convert ra/dec from catalog to ITRF coordinates
	    srcpos = radectoITRF(ateam_ra, ateam_dec, epoch, tobs_jd);
	    up = srcpos * normal > 0.1; % Sources visible at this time.
        parm.l0 = l;
        parm.m0 = m;
        parm.flux = ateam_fl;
        parm.srcname = ateam_name;
        
	    % A = exp(-(2 * pi * 1i * obs.freq / 299792458)*(posITRF* srcpos(up).'));
	    % RAteam = A * diag(ateam_fl(sel)) * A';
    else
        for ind = 1:length (parm.l0)
            parm.srcname {ind} = sprintf ('Src_%d', ind);
        end;
        up = ones (1, length (parm.l0));
    end;

    % ---------------------------------------------------------------%
    % ----------------------- Create Array --------------------------%
    % ---------------------------------------------------------------%
    % Need these parameters only for theoretical arrays
    if (strcmp (lower(parm.arrayconfig), 'lba_outer')    == 0 && ...
        strcmp (lower(parm.arrayconfig), 'lba_inner')    == 0 && ...
        strcmp (lower(parm.arrayconfig), 'lba_outer_12') == 0)
	    x = [-parm.arrayrad:parm.elemspace:parm.arrayrad]; 
        y = [-parm.arrayrad*parm.stretch:...
                            parm.elemspace:parm.arrayrad*parm.stretch];
        % Slope is always along the X dimension only.
        z = linspace (0, parm.zslope, length (x));

        % NOTE: The above array is currently at the south pole, rotate it 
        % to the location of CS002, which is where we want to place 
        % all our synthetic arrays. Rotation matrix taken from the 
        % CS002 antenna configuration file.
        cs002_rotmat = [ -0.1195950000,  -0.7919540000,   0.5987530000;
                          0.9928230000,  -0.0954190000,   0.0720990000;
                          0.0000330000,   0.6030780000,   0.7976820000];
        
        rot = [x; y; z;]' * cs002_rotmat;
        arraysampling_x = rot(:,1);
        arraysampling_y = rot(:,2);
        arraysampling_z = rot(:,3);
    end;

    fprintf (2, '<-- Simulation for timestamp %s, freq %f.\n', ...
            datestr(out.tobs), out.freq);

    
    
	deb = parm.deb;
    assert (length(parm.l0) == length (parm.m0));
    assert (length(parm.flux) == length (parm.l0));
    fprintf (1, '<-- Simulated source properties (l,m, ra, dec, Jy(frac.power),Name): \n');
    for ind = 1:length(parm.l0)
        if (parm.ateam == 1)
            fprintf (1, '   [%5.2f, %5.2f, %5.2f, %5.2f, %8.2f (%6.4f), %s]\n',...
                    parm.l0(ind), parm.m0(ind), ateam_ra(ind), ateam_dec(ind),...
                    parm.flux(ind), parm.flux(ind)/sum(parm.flux),...
                    parm.srcname{ind});
        else
            fprintf (1, '   [%5.2f, %5.2f, %8.2f (%6.4f)]\n',...
                    parm.l0(ind), parm.m0(ind), parm.flux(ind),...
                    parm.flux(ind)/sum(parm.flux));
        end;
    end 
    fprintf (1,'\n');
	l0 = parm.l0(up);
	m0 = parm.m0(up);
    flux = parm.flux(up);

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
            zpos = meshgrid (arraysampling_z, zeros (size (arraysampling_x)))';
            xpos = xpos (:);
            ypos = ypos (:);
            zpos = zpos (:);

		case 'disk'
            % Example parameter structure:
            % parm.arrayconfig='disk'; parm.l0=0.1; parm.m0=0.25;
            % parm.arrayrad=100; parm.Nelem2rad=5; parm.deb=1;

			% Create a disk of antennas, with radius being the extent of array sampling.
            [xpos, ypos] = meshgrid (arraysampling_x, arraysampling_y);
			discsel = (sqrt(xpos.^2 + ypos.^2./parm.stretch) < parm.arrayrad);
			xpos = xpos (discsel);
			ypos = ypos (discsel);
            zpos = zeros (length (xpos), 1);

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
            zpos = zeros (length (xpos), 1);

		case 'lba_outer'
			% load ('poslocal_outer_cs07w1.mat', 'poslocal');
			% load ('poslocal_outer_cs02w0.5.mat', 'poslocal');
			load ('poslocal_outer.mat', 'posITRF');
            arraysampling_x = posITRF(:,1);
            arraysampling_y = posITRF(:,2);
            remants = setdiff([1:size(posITRF,1)], parm.flagant);
			xpos = posITRF(:,1); % - posITRF(1,1);
			ypos = posITRF(:,2); % - posITRF(1,2);
            zpos = posITRF(:,3);

		case 'lba_inner'
			load ('poslocal_inner.mat', 'posITRF');
			xpos = posITRF(:,1) - posITRF(1,1);
			ypos = posITRF(:,2) - posITRF(1,2);
            zpos = posITRF(:,3);

		case 'lba_outer_12'
            load ('poslocal_afaac12_outer_0w.mat', 'poslocal');
			% load ('poslocal_afaac12_outer.mat', 'poslocal');
			xpos = poslocal(:,1); % - poslocal(1,1);
			ypos = poslocal(:,2); % - poslocal(1,2);
            zpos = poslocal(:,3);

		otherwise
			error (sprintf ('### Config %s is unimplemented!\n', ...
                    parm.arrayconfig));
	end;

    

    % Rotate the array configuration
    tmp = rotmat * [xpos(:) ypos(:)]';
    xpos = reshape (tmp (1,:), size (xpos));
    ypos = reshape (tmp (2,:), size (ypos));

	Nelem = length (xpos(:));

	if (deb > 0)
		% Show the array layout. A different colored dot for each row of 
        % elements in the array.
        figure();
		plot3 (xpos(:), ypos(:), zpos(:), '.'); 
		xlabel ('xpos (m)'); ylabel ('ypos(m)'); zlabel ('zpos(m)');
		title (sprintf ('Array configuration %s.\n', parm.arrayconfig));
	end;

    uloc = (meshgrid (xpos(:)) - meshgrid (xpos(:)).'); % U in m 
    vloc = (meshgrid (ypos(:)) - meshgrid (ypos(:)).'); % V in m
    wloc = (meshgrid (zpos(:)) - meshgrid (zpos(:)).'); % W in m
    uvdist = sqrt (uloc(:).^2 + vloc(:).^2 + wloc(:).^2);
    
    % ---------------------------------------------------------------%
    % ----------------------- Create visibilities--------------------%
    % ---------------------------------------------------------------%
    % Visibilities with an extra w-term, meant for phase-tracking
    % interferometers
    % V = sum(repmat ((parm.flux), size(uloc(:)),1) .* ...
    %       exp (-(2*pi*1i*parm.freq/299792458)*(uloc(:)*l0 + ...
    %           vloc(:)*m0 + wloc(:)*(sqrt(1-l0.^2-m0.^2) - 1))), 2);

    % Vis. as seen by arrays including a w-term
    % V = sum(repmat ((parm.flux), size(uloc(:)),1) .* ...
    %       exp (-(2*pi*1i*parm.freq/299792458)*(uloc(:)*l0 + ...
    %           vloc(:)*m0 + wloc(:)*(sqrt(1-l0.^2-m0.^2)))), 2);

    V = zeros (length(xpos));
    C = 299792458; % m/s
    posITRF = [xpos ypos zpos];
    if (isempty (parm.gain))
        parm.gain = ones (size (uloc, 1)) + 1i*zeros (size (uloc, 1));
        parm.sigman = zeros (size (uloc));
    end;

    if (isempty (parm.lofarskymodel))
        % Vis. from arrays without a w-term.
        V = sum(repmat ((flux), size(uloc(:)),1) .* ...
            exp (-(2*pi*1i*parm.freq/299792458)*(uloc(:)*l0 + vloc(:)*m0 )), 2);
        V = conj (reshape (V, [length(xpos(:)), length(ypos(:))]));
    else
        % Accumulate the visibilities due to every model source.
        for patch = 1:1 %length (rarng) % For every patch in the model.
            pos  = radectoITRF (rarng{patch}, decrng{patch}, ...
                true, tobs_jd);
            simsky_up = pos * normal > 0;
            pos (1:10,:)
            if (sum (simsky_up) == 0)
                fprintf (2, ...
                 'simvis: %s patch not visible for time %s.\n', ...
                    mod(patch).name{1}, datestr (out.tobs));
                continue;
            end;
    
    	    simsky_A = exp(-(2*pi*1i*parm.freq/C)*(posITRF*pos.'));
    	    V = V +  simsky_A*modimg{patch}*simsky_A';
        end;
    end;



    % Add noise to the visibilities
    if (parm.snr ~= 0)
        V = V + (mean(real(V(:)))/parm.snr) * randn (size (V)) + ...
                                i*(mean(imag(V(:)))/parm.snr) * randn (size (V));
    end;

    % Add the supplied antenna gains and colored noise.
    % V = V .* (parm.gain*parm.gain')/size (uloc, 1) + parm.sigman;
    
    % Another approach of first generating phasors in position coordinates.
    % Generate phasor due to location of each element
    % we = repmat ((parm.flux), size (xpos(:)), 1) .* ...
    %       (exp (-(2*pi*1i*parm.freq/299792458) *(xpos(:)*l0 + ypos(:)*m0)...
    %        + zpos(:)*sqrt(1-l0.^2-m0.^2))); 
	% V = we * we'; % Generate the visibilities for the system at hand.

    % ---------------------------------------------------------------%
    % ----------------------- Create images -------------------------%
    % ---------------------------------------------------------------%
    if (parm.fft == 0)

	    % Create an image using acm2skymap:
	    out.img_l = [-1:0.01:1];
    	out.img_m = out.img_l;
        mask = meshgrid(out.img_l).^2 + meshgrid(out.img_m).'.^2 < 1;
        out.map = acm2skyimage (V, xpos(:), ypos(:), zpos(:), out.freq, ...
                                out.img_l, out.img_m);
        out.map = flipud(rot90 (out.map)); % To match the locations of the
                                           % simulated point sources.
        scalefac = length (V(:));
    else

	    % Create a map using FFT imaging.
		gparm.type = 'pillbox';
	    gparm.lim  = 0;
	    gparm.duv = 0.5;			% Default, reassigned from freq. of obs. to
									% image just the full Fov (-1<l<1)
	    gparm.Nuv = 1500;			% size of gridded visibility matrix
	    gparm.uvpad = 1536;			% specifies if any padding needs to be added
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
    out.uloc = uloc;
    out.vloc = vloc;
    out.wloc = wloc;

    % ---------------------------------------------------------------%
    % ----------------------- Create plots if reqd-------------------%
    % ---------------------------------------------------------------%
	if (deb > 0)
        figure();
        mask = meshgrid(out.img_l).^2 + meshgrid(out.img_m).'.^2 < 1;
		% imagesc (out.img_l, out.img_m, 10*log10(real(out.map)).*mask); 
        % colorbar; axis tight;
		imagesc (out.img_l, out.img_m, (real(out.map)).*mask); 
	    axis equal
   		axis tight
   		set (gca, 'YDir', 'Normal'); % To match orientation with station images
   		set (gca, 'XDir', 'Reverse'); % To match orientation with station images

    	ylabel('South $\leftarrow$ m $\rightarrow$ North', 'interpreter', 'latex');
    	xlabel('East $\leftarrow$ l $\rightarrow$ West', 'interpreter', 'latex');
        colorbar; axis tight;
		title (sprintf ...
                ('Simulated map for %s array, max bline %.2f lambda, Ampl', ...
                parm.arrayconfig, max(uvdist(:))/(299792458/parm.freq)));
        
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
