% Script to quickly plot various parameters of visibilities.
% pep/22Nov12
% Arguments:
%	fname : Name of bin file containing visibilities to examine.
%	offset: Record offset from which to examine visibilities.
%	posfname: Antenna configuration specification for data in fname.
%	flag  : Turn on automatic flagging or not.
%   kbhit : Turn on pausing at every frame

function examinevis (fname, offset, posfname, flag, kbhit)
	[ntimes, tmin, tmax, dt] = getnrecs (fname);
	disp (sprintf ('Found %d timeslices at %f sec. resolution.\n', ntimes, dt));
	fid = fopen (fname, 'rb');
	load (posfname, 'poslocal');
    normal  = [0.598753, 0.072099, 0.797682].'; % Normal to CS002

	winsize = 10;       % NOTE: things shouldn't change over the windowsize!
	medianthresh = 2.5; % Reject timeslices with median > medianthresh*movmed;
	visamphithresh= 1.5;% Reject visibilities with median >visampthresh*median.
	visamplothresh= 0.5;% Reject visibilities with median >visampthresh*median.
	window = zeros (1, winsize);
	bline_sel = triu(ones(288)) - eye(288); % Selects upper triangle from ACM

    uloc = meshgrid (poslocal(:,1)) - meshgrid (poslocal (:,1)).';
    vloc = meshgrid (poslocal(:,2)) - meshgrid (poslocal (:,2)).';
    wloc = meshgrid (poslocal(:,3)) - meshgrid (poslocal (:,3)).';
	uvw = [uloc(:), vloc(:), wloc(:)];
	uvdist = sqrt(sum(uvw.^2, 2) - (uvw * normal).^2);
	
	% Figure management
	visimg = figure;
	visplt = figure;
	snrplt = figure;
	set(0,'Units','pixels') 
	scnsize = get(0,'ScreenSize');

	position = get(visimg,'Position');
	outerpos = get(visimg,'OuterPosition');
	borders = outerpos - position;
	edge = -borders(1)/2;
	% pos = [left bottom width height]
	pos1 = [edge, scnsize(4) * (1/2), scnsize(3)/2 - edge, scnsize(4)/2];
	set(visimg,'OuterPosition',pos1);
	pos1 = [edge+scnsize(3)/2, scnsize(4) * (1/2), scnsize(3)/2 - edge, ... 
			scnsize(4)/2];
	set(visplt,'OuterPosition',pos1);
	pos1 = [edge, scnsize(4), scnsize(3)/2 - edge, ... 
			scnsize(4)/2];
	set(snrplt,'OuterPosition',pos1);

	visnorm = zeros (1, ntimes);
	
	% Move to desired file offset.
	[acc, tobs_first, freq] = readms2float (fid, offset, -1, 288);

	% Create windows for SNR estimation on visibilities
	re_vis_win = zeros (winsize, sum(sum(bline_sel)));
	im_vis_win = re_vis_win;

	% Generate first plots to keep the axis constant.
	ind = 1;
	[acc, tobs, freq] = readms2float (fid, -1, -1, 288);
	if (isempty (tobs) == 1)
		disp ('Examinevis: End of file reached! Quitting!');
		return;
	end;

	% Flagging related.
	[uvflag, flagant] = flagdeadcorr (acc, tobs, freq, visamphithresh, ...
										visamplothresh);
	acc = acc - diag(diag(acc));
	acc (isnan(acc)) = 0; % NOTE: Norm does not work with NaN values!
	if (flag == 1)
		acc = acc .* (1-uvflag); 
	end
	visnorm (ind) = norm (acc, 'fro');
	figure (visimg);
	imagesc (10*log10(abs(acc - diag(diag(acc)))));
	% imagesc (angle(acc - diag(diag(acc))));
	title (sprintf ('Rec: %04d, Sec: %.2f (%.2f Hz)', ind+offset, tobs, freq));
	colorbar;
	disp (sprintf ('Rec:%04d, Sec: %.2f', ind+offset, tobs)); 
	
	figure (visplt);
	title (sprintf ('Rec: %04d, Sec: %.2f (%.2f Hz)', ind+offset, tobs, freq));
	subplot (1,3,1);
	plot (real (acc(:)), imag(acc(:)), '.');
	xlabel ('Real'); ylabel ('Imag');
	v  = axis; 

	subplot (1,3,2);
	plot (uvdist(:), abs (acc(:)), '.');
	xlabel ('UV Dist. (m)'); ylabel ('Visibility amp');
	uv = axis;
	
	subplot (1,3,3);
	plot (ind, visnorm(ind), 'ro');
	xlabel (sprintf ('Time (offset from %f)', tobs_first)); 
	ylabel ('Frobenius norm');
	hold on;

	% flagant = unique ([flagant missant]);

	% Fill window for statistics
	for ind = 1:winsize
		[acc, t_obs, freq] = readms2float (fid, -1, -1, 288);
		acc (isnan(acc)) = 0; % NOTE NOTE! Norm fails if members are NaNs!
		window(ind) = norm (acc - diag(diag(acc)), 'fro');
		acc_triu = acc (bline_sel(:) == 1);
		im_vis_win (ind, :) = imag (acc_triu(:));
		re_vis_win (ind, :) = real (acc_triu(:));
		normts(ind) = window(ind);
	end;
	movmed = median (window);
	re_snr = (mean(re_vis_win, 1)./std(re_vis_win));
	im_snr = (mean(im_vis_win, 1)./std(im_vis_win));


	figure (snrplt);
	subplot (1,2,1);
	plot (re_snr, 'ro');
	xlabel ('bline');
	ylabel ('Real component SNR');
	% snr = axis;

	subplot (1,2,2);
	plot (im_snr, 'bo');
	xlabel ('bline');
	ylabel ('Imag component SNR');
	% axis (snr);

	for ind = 2:ntimes
		[acc, tobs, freq] = readms2float (fid, -1, -1, 288);
		if (isempty (tobs) == 1)
			disp ('Examinevis: End of file reached! Quitting!');
			break;
		end;
		acc = acc - diag(diag(acc));
		acc (isnan(acc)) = 0; % NOTE: Norm does not work with NaN values!

		visnorm (ind) = norm (acc, 'fro');
		if (flag == 1)
			% Check for bad timeslices.
			if (visnorm(ind) > medianthresh*movmed)
				fprintf (2, '<--Discarding rec: %03d, Time: %.2f, Excess: %.2f, norm: %f, mednorm: %f\n', ...
			ind+offset, tobs, visnorm(ind)/(medianthresh*movmed), visnorm(ind), movmed);
				continue;
			end;
			acc = acc .* (1-uvflag);
		end;
			% NOTE: Ignoring bad norm before the 'continue' as well!
			window (mod(ind,winsize) + 1) = visnorm (ind); 
			acc_triu = acc (bline_sel(:) == 1);
			im_vis_win (mod(ind, winsize)+1, :) = imag (acc_triu(:)); 
			re_vis_win (mod(ind, winsize)+1, :) = real (acc_triu(:)); 
			movmed = median (window);
			re_snr = (mean(re_vis_win, 1)./std(re_vis_win));
			im_snr = (mean(im_vis_win, 1)./std(im_vis_win));
	
		figure (visimg);

		imagesc (10*log10(abs(acc - diag(diag(acc)))));
		% imagesc (angle(acc - diag(diag(acc))));
		title (sprintf ('Rec: %04d, Sec: %.2f (%.2f Hz)', ind+offset, tobs, freq));
		colorbar;
		fprintf (1, 'Rec:%04d, Sec: %.2f, norm: %f, movmed: %f.\n', ... 
				 ind+offset, tobs, visnorm(ind), movmed); 
		
		figure (visplt);
		title (sprintf ('Rec: %04d, Sec: %.2f (%.2f Hz)', ind+offset, tobs, freq));
		subplot (1,3,1);
		plot (real (acc(:)), imag(acc(:)), '.');
		axis (v);
		xlabel ('Real'); ylabel ('Imag');

		subplot (1,3,2);
		plot (uvdist(:), abs (acc(:)), '.');
		axis (uv);
		xlabel ('UV Dist. (m)'); ylabel ('Visibility amp');
		
		subplot (1,3,3);
		plot (ind, visnorm(ind), 'ro');
		plot (ind, movmed, 'bo');
		xlabel (sprintf ('Time (offset from %f)', tobs_first)); 
		ylabel ('Frobenius norm');
		hold on;

%		figure (snrplt);
%		subplot (1,2,1);
%		plot (re_snr, 'rd');
%		% hold on;
%		% axis (snr);
%		xlabel ('bline');
%		ylabel ('Real component SNR');
%
%		subplot (1,2,2);
%		plot (im_snr, 'bd');
%		% hold on;
%		% axis (snr);
%		xlabel ('bline');
%		ylabel ('Imag component SNR');
		if (kbhit == 1)
			pause; %  (0.01);
		end;
	end
	% figure;
	% plot (visnorm);

