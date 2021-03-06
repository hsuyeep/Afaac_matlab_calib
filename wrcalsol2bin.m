% Script to write out Calibration solutions as floats to binary file.
% NOTE: Always writes out 288 sized solutions,regardless of flagged antennas.
% The unused rows and columns are filled with 0s, based on the antenna 
% number of the dipoles which are flagged, specified in flagants.
% pep/18Oct12
% Arguments:
%	fid : File id of solution file to write to.
%   sol : Solution structure as generated by pelican_sunAteamsub.m etc.

function wrcalsol2bin (fid, sol)
	Nelem = 288; 

	if (fid < 0 || isempty(sol))
		err = MException('wrcalsol2bin:fid < 0', ...
        		'file id or sol. inconsistent!');	
		throw (err);
	end;
	
	% Write out a single record, corresponding to a single time/freq bin.
	% Assumes flagants are appropriately identifiable in the input data.
	% NOTE: Should move gainsol size consolidation with flagged antennas to 
	% layer above at some point... 
	gainsol = complex (zeros (1, Nelem), zeros (1, Nelem));
	% gainsol(gainmask == 0) = sol.gainsol; % NOTE: Gains always have
	% rodata.Nelem elements, with flagged antennas having 0 gain.
	gainsol = sol.gainsol;
	count = fwrite (fid, sol.tobs, 'double');
	count = fwrite (fid, sol.freq, 'double');

	% The flagged antennas
	count = fwrite (fid, length (sol.flagant), 'float32');
	count = fwrite (fid, sol.flagant, 'float32');

	% NOTE:store len of gainsol array, ideally should always be Nelem.
	count = fwrite (fid, single(length (gainsol)), 'float32'); 
	count = fwrite (fid, single(real(gainsol)), 'float32');
	count = fwrite (fid, single(imag(gainsol)), 'float32');

	% Statistics of cal_ext_stefcal effort in estimation
	count = fwrite (fid, single (sol.calext_iters), 'float32');
	count = fwrite (fid, single (sol.pinv_sol), 'float32');

	% Statistics of gainsolv effort in estimation
	count = fwrite (fid, single (sol.stefcal_iters), 'float32');
	count = fwrite (fid, single (sol.stefcal_tol), 'float32');

	% NOTE: store len of model sky sources whose flux has been extracted.
	count = fwrite (fid, single(length (sol.sigmas)), 'float32'); 
	count = fwrite (fid, single(sol.sigmas), 'float32');

	% NOTE: store len of coordinates of extracted sources. Currently they can be
	% different from the number of srcs. whose flux has been extracted (above)
	count = fwrite (fid, single(length (sol.thsrc_cat)), 'float32'); 
	count = fwrite (fid, single(sol.thsrc_cat), 'float32');
	count = fwrite (fid, single(sol.phisrc_cat), 'float32');
	count = fwrite (fid, single(sol.thsrc_wsf), 'float32');
	count = fwrite (fid, single(sol.phisrc_wsf), 'float32');

	% NOTE: Store length of sigman, which only includes unflagged antennas!
	count = fwrite (fid, single (length (sol.sigman(:))), 'float32'); 
	count = fwrite (fid, single (real(sol.sigman(:))), 'float32');
	count = fwrite (fid, single (imag(sol.sigman(:))), 'float32');
