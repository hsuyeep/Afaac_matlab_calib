% Script to generate statistics of observed visibilities available in a .bin
% file as generated by ms2float.py
% pep/18Oct12

function [sig, mu, histre, histim] = genvisstats (fname, ntimes)
	if ntimes == -1
		[ntimes, tmin, tmax, dt] = getnrecs (fname);
	end
	fid = fopen (fname, 'rb');
	bins = 16;
	binsby2p1 = bins/2 + 1;
	histmin = -0.05; histmax = +0.05; % chosen based on manual check of data.
	histdt = (histmax-histmin)/bins;
	[acc, tobs, freq] = readms2float (fid, -1, -1, 288);
	au = triu (acc) - diag(diag(acc));
	au_vec = au (au ~= 0);
	fseek (fid, 0, 'bof');
	nelem = size (acc, 2);
	nbl = nelem * (nelem - 1)/2;

	% Initialize accumulators
	bl = au_vec; % Store for sum(visibilities);
	bl2 = (real(au_vec).^2 + i*imag(au_vec).^2); % Acc squares.
	histre = zeros (bins, length(au_vec));
	histim = zeros (bins, length(au_vec));

	for ts = 1:ntimes
		[acc, tobs, freq] = readms2float (fid, -1, -1, 288);
		au = triu (acc) - diag(diag (acc)); % Eliminate autocorrelations.
		au_vec = au (au ~= 0);
		bl = bl + au_vec;     				% Accumulate complex visibilities
		bl2 = bl2 + (real(au_vec).^2 + i*imag(au_vec).^2); % Acc squares.
		indre = round(real(au_vec) ./ histdt) + binsby2p1;
		indre (indre < 1) = 1; indre(indre > bins) = bins;
		indim = round (imag(au_vec) ./ histdt) + binsby2p1;
		indim (indim < 1) = 1; indim(indim > bins) = bins;
		for b = 1:nbl;
			histre (indre(b), b) = histre (indre(b), b) + 1;
			histim (indim(b), b) = histim (indim(b), b) + 1;
		end;
	end;		

	% Compute stats.
	mu = bl ./ ntimes;
	% var = Mean of the squares - sq. of the means
	sig = (bl2 ./ ntimes) - (real(mu).^2 + i*imag(mu).^2); 
