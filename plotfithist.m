% Function to overplot a fitted gaussian over the histogram of data.
% Slightly modified from the matlab online help section for data analysis.
% pep/26Jul13
% Arguments:
%	vec: Vector of data
%  nbin: Number of bins required in the histogram
%   hdl: Plot handle, usually obtained by calling 'gca'

function redchisq = plotfithist (vec, nbin, hdl)
	if (isempty (vec)) 
		fprintf (2, 'plotfithist: Empty vector!\n'); return;
	end;
	[bin_cnt, bin_loc] = hist (vec, nbin);
	binwid = bin_loc(2) - bin_loc (1);
	histarea = binwid*(sum(bin_cnt));
	hist (vec, nbin);
	hold on;
	m = mean (vec(isnan(vec) ~= 1)); s = std (vec(isnan(vec) ~= 1));
	pdf = histarea* (1/(s*sqrt(2*pi)))*exp (-(bin_loc - m).^2/(2*s.^2));
	plot (bin_loc, pdf, 'r', 'LineWidth', 2);
	hold off;
	xlabel ('bin'); ylabel ('Freq'); title ('Hist. with gauss. fit');
	text (bin_loc(int32(3*nbin/4)), max(bin_cnt), ...
		  sprintf('mean: %.4f\nstd : %.4f', m, s));

	% Compute chisq. to determine how good a fit the model gaussian is to 
	% the generated histogram. See Sec. 4.3, Bevington
	sel = (isnan (pdf) == 0) & (bin_cnt > 0);
	chisq = sum ((bin_cnt(sel) - pdf(sel)).^2./bin_cnt(sel));
	redchisq = chisq/(nbin-1);  % Reduced chisq, with nbin-1 deg. of freedom.

