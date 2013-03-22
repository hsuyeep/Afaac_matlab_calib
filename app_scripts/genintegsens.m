% Process to generate a graph of increment of sensitivity with integration
% time, after rephasing over a specifiable window of records. 
% pep/25Nov12
	% clear all; close all;
	% fname = '~/WORK/AARTFAAC/rawdat/SB000_ch30-35_5sec_3hr_cal.bin';
	% integ_secs = 5;
	debug = 3;
	nfacet = 3, facetsize = 256;

	fname = '~/WORK/AARTFAAC/Reobs/11Jul12/LBA_OUTER_BAND60/SB001_LBA_OUTER_4min_1ch_cal.bin';
	integ_secs = 1;
    [ntslices, tmin, tmax, dt] = getnrecs (fname);
	% ntslices = 250;
	disp (sprintf ('Found %d records at %f sec. resolution.', ntslices, dt));
	fid = fopen (fname, 'rb');
	[acct0, tobs, freq] = readms2float (fid, -1, -1);
	tint = 3:8:119;
	acc = zeros ([size(acct0), tint(end)]);
	tobs = zeros (1, tint(end));
	load ('poslocal.mat', 'poslocal','posITRF');
	load ('srclist3CR.mat');
	uloc = meshgrid (poslocal(:,1)) - meshgrid (poslocal(:,1)).';
	vloc = meshgrid (poslocal(:,2)) - meshgrid (poslocal(:,2)).';
	% imaging parameters
	duv = 2.5; Nuv = 700; uvpad = 768; % Each facet is 512x512 pixels
	mask = [];
	dr = zeros (int32(ntslices/tint(end)), length(tint));
	sig = zeros (int32(ntslices/tint(end)), length(tint));
	dr_diff = zeros (int32(ntslices/tint(end)), length(tint));
	sig_diff = zeros (int32(ntslices/tint(end)), length(tint));
	% sig = [];
	dr_ind = 1;

	plthdl = figure;
	imghdl = figure;
	diffhdl = figure;

	% for trec = 1:tint(end):ntslices
	for trec = 1:tint(end):tint(end)
		dr_ind = 1;

		% NOTE: This loop reads in ACMs for the largest window defined, and 
		% further code works on subsets of the full window.
		for ind = 1:tint(end)
				[acct0, t, freq] = readms2float (fid, -1, -1);
				if (isempty (t)) 
					break;
				end;
				% Skipping bad data, for LBA_OUTER_BAND60/*4min*cal.bin
				badind = trec + ind ;
				if (badind==35 || badind == 68 || badind == 80 || badind == 151 || badind == 234)
					disp(sprintf ('Skipping timeslice %d, time %f.', ind, t));
					acc(:,:,ind) = acc(:,:,ind-1);
					tobs(ind) = t;
					continue;
				end
				acc(:,:,ind) = acct0;
				tobs(ind) = t;
				disp (['Time: ' num2str(tobs(ind))]);
		end;
		if (isempty (t))
			disp ('End of file reached');
			break;
		end;
			
		% Loop over different integration time windows.
		disp(['Operating on time range (' num2str([tobs(1), tobs(end)]) ')']);
		for jind = 1:length(tint);
	
			% Integrate visibilities over the specified window (tint(jind)),
			% Rephasing is done to a time corresponding to the center of the 
			% tobs array. The integrated visibilities are returned in 
			% integ_acc.
			[reph_acc, integ_acc, integ_tobs] = integvis (acc, tint(jind), ...
												tobs, freq, posITRF, 0);
		
			%-------- Derive statistics from integrated visibilities -------
			% Choose imager for integrated visibilities.
			[radecmap1, calmap, calvis, l, m] = ... 
			    fft_imager_sjw_radec(integ_acc (:),...
		   	       uloc(:), vloc(:), duv, Nuv, uvpad, integ_tobs, freq, 0);

			% Image current timeslice. Generate a mosaic image.
	   		% [calmap, l, m] = genmosaic(integ_acc, integ_tobs, freq, ...
		 %					 nfacet, facetsize, uloc(:), vloc(:), posITRF, 0);
			if isempty (mask)
		    	% mask = NaN (length (l));
		    	mask = zeros (length (l));
				mask ((meshgrid (l).^2 + (meshgrid(l).').^2) < 1) = 1;
			end;
			map = abs (calmap) .* mask;
			map_nonan = abs(calmap);
			map_nonan (isnan(map_nonan)) = 0;
			[zen_dr, zen_sig] = getimagedr (map, 64, 3);
			dr (trec, dr_ind) = zen_dr;
			sig (trec, dr_ind) = zen_sig;

			% To plot map made from integrated visibilities
			if debug > 2
				figure (imghdl);
				imagesc(l,m,map);
				title (sprintf  ('Integ: %d, dr:%f, sig:%f', ... 
						tint(jind), zen_dr, zen_sig));
%				figure (plthdl);
%				hist (map_nonan, 100);
				% pause
				% fseek (fid, 0, -1);
			end

			%-------- Statistics from rephased image diffs over window -------
			tmpacc = reph_acc (:, :, 1); % First rephased ACM of this window.
			[radecmap1, calmap, calvis, l, m] = ... 
			    fft_imager_sjw_radec(tmpacc (:),...
		   	       uloc(:), vloc(:), duv, Nuv, uvpad, tobs(1), freq, 0);
			map = abs (calmap) .* mask;
			prevmap = map;
			diffmap = zeros (size (prevmap));

			% Operate on each rephased ACM of this window.
			for ind = 2:tint(jind)
				tmpacc = reph_acc (:, :, ind); 
				[radecmap1, calmap, calvis, l, m] = ... 
			   	 fft_imager_sjw_radec(tmpacc (:),...
		   	   	    uloc(:), vloc(:), duv, Nuv, uvpad, tobs(ind), freq, 0);
				map = abs (calmap) .* mask;
				diffmap = diffmap + (map - prevmap);
				prevmap = map;
				% To plot maps made from instantaneous, rephased ACMs.
%				if debug > 2
%					figure (imghdl);
%					imagesc(l,m,map);
%					title (sprintf  ('Instantaneous map: tobs = %f', ... 
%							tobs(ind)));
%					figure (diffhdl);
%					imagesc(l,m,diffmap);
%					title (sprintf  ('Snapshot of cumulative difference map: tobs = %f', ... 
%							tobs(ind)));
%				end
			end
			diffmap = diffmap ./tint(jind);
			[zen_dr, zen_sig] = getimagedr (diffmap, 64, 3);
			dr_diff (trec, dr_ind) = zen_dr;
			sig_diff (trec, dr_ind) = zen_sig;
			% Histogram of averaged difference images over window.
			% figure (plthdl);
			% hist (diffmap(:), 100);
			title (sprintf ('Histogram of average difference image over window size %d', tint(jind)));


			dr_ind = dr_ind + 1;
			imgmax = max(max(map));
			imgmin = min(min(map));     

		end	
		figure; %  (plthdl);
		subplot (1,3,1); 
		plot (integ_secs*tint, dr(trec,:), 'b.-'); 
		title ('Dynamic range Vs. integration time');
		xlabel ('Integration time (secs)');
		ylabel ('Dynamic range (linear)');

		subplot (1,3,2); 
		% plot (integ_secs*tint, sig(trec,:), 'b.-'); 
		plot (integ_secs*tint, sig(trec,:), 'b.-'); 
		title ('rms noise Vs. integration time');
		xlabel ('Integration time (secs)');
		ylabel ('sigma');

		subplot (1,3,3); 
		% plot (integ_secs*tint, sig(trec,:), 'b.-'); 
		plot (integ_secs*tint, sig_diff(trec,:), 'b.-'); 
		title ('rms noise Vs. integration time');
		xlabel ('Integration time (secs)');
		ylabel ('sigma');
		pause;
	end
	figure (plthdl);
	subplot (1,2,1); 
	plot (integ_secs*tint, mean(dr), 'b.-'); 
	title ('Dynamic range Vs. integration time');
	xlabel ('Integration time (secs)');
	ylabel ('Dynamic range (linear)');

	subplot (1,2,2); 
	% plot (integ_secs*tint, sig(trec,:), 'b.-'); 
	plot (integ_secs*tint, mean(sig), 'b.-'); 
	title ('rms noise Vs. integration time');
	xlabel ('Integration time (secs)');
	ylabel ('sigma');
	pause;

