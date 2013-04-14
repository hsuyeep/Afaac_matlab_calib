% Script to generate statistics on the distribution of baselines for various
% AARTFAAC configurations.
% pep/09Apr13

function blinedist (array_config)
	% Array config is one of 'LBA_OUTER', '
	allowed_arr_confs = ... 
		{'LBA_OUTER', 'LBA_INNER', 'LBA_X', 'LBA_Y', 'LBA_SPARSE', ...
		 'HBA_DUAL'};
	C = 299792458;					% Speed of light, m/s.
	radec = 0;
    duv = 2.5;						% Default, reassigned from freq. of obs. to
									% image just the full Fov (-1<l<1)
    Nuv = 500; %1000                % size of gridded visibility matrix
    uvpad = 512; %1024              % specifies if any padding needs to be added
	Nelem = 288;
	acc = ones (288);
	acc = acc - diag (diag (acc));

	% Load various antenna configurations
	switch array_config
		case 'LBA_OUTER'
			load ('poslocal_outer.mat', 'poslocal');
			freq = 60000000;
			nameind = 1;

		case 'LBA_INNER'
			load ('poslocal_inner.mat', 'poslocal');
			freq = 60000000;
			nameind = 2;
		
		otherwise
			fprintf (2, 'blinedist: Unknown array configuration!');

	end;
	

	% Generate baselines at different freq (hence spatial frequencies)
	uloc = meshgrid (poslocal (:,1)) - meshgrid (poslocal (:,1)).';
	vloc = meshgrid (poslocal (:,2)) - meshgrid (poslocal (:,2)).';
	wloc = meshgrid (poslocal (:,3)) - meshgrid (poslocal (:,3)).';
	uvw = [uloc(:), vloc(:), wloc(:)];
    normal  = [0.598753, 0.072099, 0.797682].'; % Normal to CS002
	uvdist = sqrt(sum(uvw.^2, 2) - (uvw * normal).^2);

	% Generate distribution of baselines 
	figure;
	subplot (1,2,1);
	hist (uloc(:), 100);
	title ('Baseline u components');
	subplot (1,2,2);
	hist (vloc(:), 100);
	title ('Baseline v components');

	% Generate PSF profiles for various weighting and tapering schemes.
	% Natural weighting
	taper = acc(:);
	% Put in taper according to visibility cutoff in calib. code.
	% mask = reshape(uvdist, [Nelem, Nelem]) < ... 
	% 			min([mrestriction * (C / freq), maxrestriction]);
	
	figure;
   	[radecmap, map, calvis, l, m] = ... 
		  fft_imager_sjw_radec (acc(:), uloc(:), vloc(:), ... 
					duv, Nuv, uvpad, 0, freq, radec);
	mesh (l, m, double(20*log10(abs(map)/max(abs(map(:))))))
	title (sprintf ('PSF in dB for %s array configuration at %.2f MHz, Natural', ... 
					 char(allowed_arr_confs(nameind)), freq/1e6));

	% Uniform: weigh according to density of baselines
	cellsize = 10; % m, all baselines within this interval
	uvec = uloc (:); vvec = vloc (:);
	accvec = acc (:);
	for uind = 1:max(uvec)/cellsize
		usel = (abs(uvec) > (uind-1)*cellsize) & (abs(uvec) < uind*cellsize);
		for vind = 1:max(vvec)/cellsize
			vsel = (abs(vvec) > (vind-1)*cellsize) & (abs(vvec)<vind*cellsize);
			sel = usel & vsel;
			w = sum (sel); % Weight = total number of visibilities in range.
			accvec (usel & vsel) = 1/w;
		end;
	end;
	acc_weighted= reshape (accvec, [Nelem, Nelem]);

	% Handle tapering.
	figure;
   	[radecmap, map, calvis, l, m] = ... 
		  fft_imager_sjw_radec (acc_weighted(:), uloc(:), vloc(:), ... 
					duv, Nuv, uvpad, 0, freq, radec);
	mesh (l, m, double(20*log10(abs(map)/max(abs(map(:))))))
	title (sprintf ('PSF in dB for %s array configuration at %.2f MHz, Uniform', ... 
					 char(allowed_arr_confs(nameind)), freq/1e6));
