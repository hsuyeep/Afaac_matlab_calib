% Script to display images from multiple subbands as a matrix.
% Based on genfftimage.m
% pep/22Nov12

function multibandimg (ntslices)
	sband0 = '~/WORK/AARTFAAC/Reobs/11Jul12/LBA_OUTER_BAND_SPREAD/full_conv/SB000_LBA_OUTER_SPREAD_1ch_1_convcal.bin';
    sband1 = '~/WORK/AARTFAAC/Reobs/11Jul12/LBA_OUTER_BAND_SPREAD/full_conv/SB001_LBA_OUTER_SPREAD_1ch_1_convcal.bin';
	sband2 = '~/WORK/AARTFAAC/Reobs/11Jul12/LBA_OUTER_BAND_SPREAD/full_conv/SB002_LBA_OUTER_SPREAD_1ch_1_convcal.bin';
	sband3 = '~/WORK/AARTFAAC/Reobs/11Jul12/LBA_OUTER_BAND_SPREAD/full_conv/SB003_LBA_OUTER_SPREAD_1ch_1_convcal.bin';
	sband4 = '~/WORK/AARTFAAC/Reobs/11Jul12/LBA_OUTER_BAND_SPREAD/full_conv/SB004_LBA_OUTER_SPREAD_1ch_1_convcal.bin';

	fids = zeros (1,5);
	fids(1) = fopen (sband0, 'rb');
	fids(2) = fopen (sband1, 'rb');
	fids(3) = fopen (sband2, 'rb');
	fids(4) = fopen (sband3, 'rb');
	fids(5) = fopen (sband4, 'rb');
	
	radec = 0;
    duv = 2.5;						% Default, reassigned from freq. of obs. to
									% image just the full Fov (-1<l<1)
    Nuv = 500; %1000                % size of gridded visibility matrix
    uvpad = 512; %1024              % specifies if any padding needs to be added
	hdl = figure;

	fid = zeros (1, 5);
	for ind = 1:5
		% Obtain image parameters from visibility file.
		[acc, t_obs, fr] = readms2float (fids(ind), -1, -1, 288);
		freq(ind) = fr;
		lambda(ind) = 299792458/freq(ind); 		% in m.
		duv(ind) = lambda(ind)/2;
		% dimensionless, in dir. cos. units
    	dl(ind) = (299792458/(freq(ind) * uvpad * duv(ind))); 
	end
	
	station = 0;
	% For imaging from a single station
	% acc_tmp = acc;
	% acc = zeros (48, 48);
	% acc = acc_tmp (1:48, 1:48);
	% acc = acc_tmp (49:96, 49:96);
	% acc = acc_tmp (97:144, 97:144);
    
    % NOTE: Total imaged Field of View is determined by the visibility 
	% grid-spacing, duv.
    lmax = dl * uvpad / 2;
    % l = [-lmax:dl:lmax-1e-3];
    % m = l;  % Identical resolution and extent along m-axis

    % Local horizon based coordinates of array in ITRF
    load ('poslocal.mat', 'posITRF', 'poslocal'); 

	% load 3CR catalog
	load ('srclist3CR.mat');

    % Generate uv coordinates in local horizon coord. system, needed for imaging
    uloc = meshgrid (poslocal(:,1)) - meshgrid (poslocal (:,1)).';
    vloc = meshgrid (poslocal(:,2)) - meshgrid (poslocal (:,2)).';
	
% 	for ts = 1:ntslices
		for ind = 1:5
			[acc, t_obs, fr] = readms2float (fids(ind), -1, -1, 288);
			if (isempty (acc) == false)
				disp (['Time/Freq: ' num2str(t_obs) ' ' num2str(freq)]); 
				% Image current timeslice.
	   			[radecmap, calmap, calvis, l, m] = ... 
					fft_imager_sjw_radec (acc(:), uloc(:), vloc(:), ... 
									duv(ind), Nuv, uvpad, t_obs, fr, radec);
		
				subplot (2,3,ind)
				subimage (l,m,abs(calmap));
	       		% imagesc(l, m, calmap);
%			if isempty (caxisrng) == 0
%				caxis (caxisrng);
%			end;
%			[sel, sel_l, sel_m] = overplot3cr (t_obs, srclist3CR, 20, hdl);
%	        set(gca, 'FontSize', 16);
			
	        title(['Time: ' num2str(t_obs) ', Freq: ' num2str(freq(ind))]);
	        axis equal
	        axis tight
	        xlabel('South \leftarrow m \rightarrow North');
	        ylabel('East \leftarrow l \rightarrow West');
	        set(colorbar, 'FontSize', 16);
		  end
	
%			% Read in next visibility set.
%			[acc, t_obs, freq] = readms2float (fin, -1, -1);
%			if (isempty (acc))
%				disp ('End of file reached!');
%				break;
%			end;
%			pause (0.01);
%		else
%			disp ('File end reached!');
%		end;
	end;
