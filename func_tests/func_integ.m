% Script to integrate visibilities after rephasing them to a common phase 
% center/
% pep/15Nov12
% Arguments:
%  acc: 3D cube of visibilities, in increasing time order.
%  tobs: Corresponding observing times for each acm.
% Returns:
%  integ_acc: The integrated visibilities, after rephasing
%  integ_tobs: The time corresponding to the center of the integrated data range

% function [integ_acc, integ_tobs] = integvis (acc, tobs)
function [integ_acc, integ_tobs] = func_integ(tint)
	% close all;
	% fname = '~/WORK/AARTFAAC/rawdat/SB000_ch30-35_5sec_3hr_cal.bin';
	fname = '~/WORK/AARTFAAC/Reobs/11Jul12/LBA_OUTER_BAND60/SB001_LBA_OUTER_4min_1ch_cal.bin';
	debug = 2;
	[ntimes, tmin, tmax, dt] = getnrecs (fname);
	disp (sprintf ('Found %d recs at %f resolution.', ntimes, dt));
	fid = fopen (fname, 'rb');
	[acct0, tobs, freq] = readms2float (fid, -1, -1, 288);
	acc = zeros ([size(acct0), tint]);
	tobs = zeros (1, tint);

	load ('poslocal.mat', 'poslocal','posITRF');
	load ('srclist3CR.mat');
	uloc = meshgrid (poslocal(:,1)) - meshgrid (poslocal(:,1)).';
	vloc = meshgrid (poslocal(:,2)) - meshgrid (poslocal(:,2)).';
	
	% imaging parameters
	duv = 2.5; Nuv = 700; uvpad = 768; % Each facet is 512x512 pixels
	nfacet = 3, facetsize = 256;

	for ind = 1:tint
		[acct0, t, freq] = readms2float (fid, -1, -1, 288);
		acc(:,:,ind) = acct0;
		tobs(ind) = t;
		disp (['Time: ' num2str(tobs(ind))]);
	end;

	ntslices = size (acc, 3);
	integ_acc = zeros (size (acc (:,:, 1)));
	integ_tobs = 0;
	mosimg = figure;
	mosplt = figure;
	if debug > 2
		rephimg = figure;
		rephpro = figure;
	end
	mask = [];

for nrec = 1:ntslices:ntimes
	integ_tobs = 0;
	for ind = 1:ntslices
		disp (['Rephasing from ' num2str(tobs(ind)) ' to ' num2str(tobs(length(tobs)))]);
		reph_acc = rephasetime (acc(:,:,ind), tobs(ind), tobs(length(tobs)), ... 
								freq, posITRF);
		integ_acc = integ_acc + reph_acc;
		integ_tobs = integ_tobs + tobs(ind);
		% To image each rephased ACM, or the integrated ACM for every 
		% integration.
		if debug > 2
	    	[radecmap1, rephmap, calvis, l, m] = ... 
				fft_imager_sjw_radec(reph_acc(:), ...
					 uloc(:), vloc(:), duv, Nuv, uvpad, tobs(ind), freq, 0);
	    	if isempty (mask)
		     	% mask = NaN (length (l));
		     	mask = zeros (length (l));
	   		  	mask ((meshgrid (l).^2 + (meshgrid(l).').^2) < 1) = 1;
	   		end;
			figure(rephimg);
	    	map = abs (rephmap) .* mask;
			imagesc(l,m,map);
   	 		[zen_dr, zen_sig] = getimagedr (map, 64, 3);
		    imgmax = max(max(map));
			imgmin = min(min(map));     
			title(sprintf ('%f: DR: %3.0f, sig: %3.2f, max: %5.2f, min:%5.2f',...
				 integ_tobs, zen_dr, zen_sig, imgmax, imgmin));
			map_nonan = abs(rephmap);
			map_nonan (isnan(map_nonan)) = 0;
			figure (rephpro);
			subplot (2,2,1);
			% plot (l, sum(abs(rephmap), 1));
			plot (l, sum(map_nonan, 1));
			subplot (2,2,2);
			plot (m, sum(map_nonan, 2));
			subplot (2,2,3);
			hist (map_nonan(:), 100);
			pause(0.1);
		end
	end;	
	integ_tobs = integ_tobs/ntslices; 
	integ_acc = integ_acc / ntslices;
	
	% Image current timeslice. Generate a zenith pointing image.
    % [radecmap1, calmap, calvis, l, m] = fft_imager_sjw_radec(integ_acc (:),...
    %                 uloc(:), vloc(:), duv, Nuv, uvpad, integ_tobs, freq, 0);

	% Image current timeslice. Generate a mosaic image.
    [calmap, l, m] = genmosaic(integ_acc, integ_tobs, freq, nfacet, facetsize, uloc(:), vloc(:), posITRF, 0);
	if isempty (mask)
     	% mask = NaN (length (l));
     	mask = zeros (length (l));
	  	mask ((meshgrid (l).^2 + (meshgrid(l).').^2) < 1) = 1;
	end;
    
    map = abs (calmap) .* mask;
	map_nonan = map; %  abs(rephmap);
	map_nonan (isnan(map_nonan)) = 0;
	figure (mosimg);
	imagesc(l,m,map);
	[sel, sel_l, sel_m] = overplotcat (integ_tobs, srclist3CR, 50, mosimg, true);
	colorbar;
    [zen_dr, zen_sig] = getimagedr (map, 64, 3);
    imgmax = max(max(map));
    imgmin = min(min(map));     
	title(sprintf ('%f:(%d slices). DR: %3.0f, sig: %3.2f, max: %5.2f, min:%5.2f',...
			 integ_tobs, tint, zen_dr, zen_sig, imgmax, imgmin));
	figure (mosplt);
	subplot (2,2,1);
	% plot (l, sum(abs(rephmap), 1));
	plot (l, sum(map_nonan, 1));
	subplot (2,2,2);
	plot (m, sum(map_nonan, 2));
	subplot (2,2,3);
	hist (map_nonan(:), 100);

	for ind = 1:tint
		[acct0, t, freq] = readms2float (fid, -1, -1, 288);
		if (isempty (t) == 1)
			break;
		end;
		acc(:,:,ind) = acct0;
		tobs(ind) = t;
		disp (['Time: ' num2str(tobs(ind))]);
	end;
	if (isempty (t) == 1)
		disp ('EoF reached! Bailing out...');
		break;
	end;
	pause;
end
