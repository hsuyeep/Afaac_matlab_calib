% Script to integrate visibilities after rephasing them to a common phase 
% center, defined by the center of tobs array.
% pep/15Nov12

% Arguments:
%  acc: 3D cube of visibilities, in increasing time order.
%  tobs: Corresponding observing times for each acm.
% Returns:
%  integ_acc: The integrated visibilities, after rephasing
%  integ_tobs: The time corresponding to the center of the integrated data range

function [reph_acc, integ_acc, integ_tobs] = integvis (acc, ntslices, tobs, freq, posITRF, debug)
	% ntslices = size (acc, 3);
	integ_acc = zeros (size (acc (:,:, 1)));
	reph_acc = zeros(size(acc));
	integ_tobs = 0;
	% t_reph = tobs (int32(length(tobs)/2));
	t_reph = tobs (int32(ntslices/2));

	% Debug imaging.
	if debug == 1
		rephimg = figure;
		rephpro = figure;
		load ('poslocal.mat', 'poslocal','posITRF');
		load ('srclist3CR.mat');
		uloc = meshgrid (poslocal(:,1)) - meshgrid (poslocal(:,1)).';
		vloc = meshgrid (poslocal(:,2)) - meshgrid (poslocal(:,2)).';
		% imaging parameters
		gparm.type = 'pillbox';
		gparm.duv = 2.5; 
		gparm.Nuv = 700;
		gparm.uvpad = 768; 
		gparm.fft = 1;

	end

	mask = [];

	for ind = 1:ntslices
		disp (['Rephasing from ' num2str(tobs(ind)) ' to ' num2str(t_reph)]);
		reph_acc (:, :, ind) = rephasetime (acc(:,:,ind), tobs(ind), ... 
								t_reph, freq, posITRF);


		integ_acc = integ_acc + reph_acc(:,:,ind);
		integ_tobs = integ_tobs + tobs(ind);

		if debug == 1
			% Debug imaging
			tmpacc = reph_acc(:,:,ind);
    		[radecmap1, rephmap, calvis, l, m] = fft_imager_sjw_radec(tmpacc(:), ...
					 uloc(:), vloc(:), gparm, [], [], tobs(ind), freq, 0);
    		if isempty (mask)
	     		% mask = NaN (length (l));
	     		mask = zeros (length (l));
   		  		mask ((meshgrid (l).^2 + (meshgrid(l).').^2) < 1) = 1;
   			end;
	    	map = abs (rephmap) .* mask;
	    	[zen_dr, zen_sig] = getimagedr (map, 64, 3);
		    imgmax = max(max(map));
			imgmin = min(min(map));     
			figure(rephimg);
			imagesc(l,m,map);
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
	

    
	if debug == 1
    	[radecmap1, calmap, calvis, l, m] = fft_imager_sjw_radec(integ_acc (:),...
                    uloc(:), vloc(:), gparm, [], [], tobs, freq, 0);
	    map = abs (calmap) .* mask;
		map_nonan = abs(rephmap);
		map_nonan (isnan(map_nonan)) = 0;
		figure;
		imagesc(l,m,map);
	    [zen_dr, zen_sig] = getimagedr (map, 64, 3);
	    imgmax = max(max(map));
	    imgmin = min(min(map));     
		title(sprintf ('%f:(%d slices). DR: %3.0f, sig: %3.2f, max: %5.2f, min:%5.2f',...
				 integ_tobs, ntslices, zen_dr, zen_sig, imgmax, imgmin));
		figure ;
		subplot (2,2,1);
		% plot (l, sum(abs(rephmap), 1));
		plot (l, sum(map_nonan, 1));
		subplot (2,2,2);
		plot (m, sum(map_nonan, 2));
		subplot (2,2,3);
		hist (map_nonan(:), 100);
	end
