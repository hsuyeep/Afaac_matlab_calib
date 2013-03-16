% Script to write out ACMs as floats in the format accepted by readms2float.m
% NOTE: Always writes out 288x288 ACM,regardless of flagged antennas. The 
% correct rows and columns are filled with 0s, based on the antenna number of 
% the dipoles which are flagged, specified in flagants.
% pep/18Oct12

% Arguments:
%	fid: File id of the binary file to which calibrated visibilities should be
%		 written.
%	acc: Calibrated ACM to write to file. NOTE: Size INcludes flagged ants, which
%		 have NaN as values!
% flagants: Vector containing indices of flagged antennas.
%	tobs: The MJD sec. time corresponding to this calibrated ACM.
%	freq: The Freq. in Hz. corresponding to this calibrated ACM.
%


function wracm2bin (fid, acc, flagants, tobs, freq)
	Nelem = 288; 

%	if (size (acc, 1) + length(flagants) > Nelem)
%		disp ('wracm2bin: Inconsistency in number of elements in ACM!');
%		disp (['Flagants: ' num2str(flagants) ' size (acc):' num2str(size(acc))]);
%		return;
%	end;
    nblines = Nelem * (Nelem + 1)/2;
    recsize = 8*(2+nblines); % Bytes, assumed double =8, float=4 bytes.

	% This always gives 288*289/2 complex values.
	b = acc (triu(ones (Nelem)) == 1); 

	% NOTE: This changes the size of visibilities to write out, in case there 
	% are genuine visibilities with 0 values! Subtle bug!
	% au = triu (acc);
	% b = au (au ~= 0); 

	% NOTE: Need to interleave real and imaginary values explicitly
	b_re = real (b); b_im = imag (b);
	interleaved = zeros (1, 2*length (b_re));
    for bl = 1:length (b_re)
		interleaved (2*bl-1) = b_re(bl);
		interleaved (2*bl)   = b_im(bl);
	end;
	fprintf (1, 'wracm2bin: Writing tobs: %.2f, Freq. %.2f.\n', tobs, freq);
	count = fwrite (fid, tobs, 'double');
	count = fwrite (fid, freq, 'double');
	count = fwrite (fid, interleaved, 'float32');
