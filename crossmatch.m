% Script to cross match source across catalogs. Can be used to find a matching
% source across catalogs.
% pep/20Jul13
% Arguments:
%	cat1 : The first catalog, usually the 3CR.
%	ind1 : List of sources (as indices) to crossmatch.
%	ep1  : The epoch of cat1, true => B1950, false => J2000
%	cat2 : The second catalog, usually VLSS
%	ep2  : The epoch of cat2, same as above.
% Returns:
%	ind2 : The indices of sources in cat2 which match cat1. If a particular
%			source is not found, its index is set to 0.
%   cat1a2k, cat1d2k: The J2000 coordinates of cat1. sources.

function [ind2, cat1a2k, cat1d2k] = crossmatch (cat1, ind1, ep1, cat2, ep2, debug)

	% Conversion from B1950 to J2000.
	precMat = precessionMatrix(JulianDay(datenum(1950, 1, 1, 0, 0, 0)));

	% Acceptable error radius, in rad.
	errr = 0.0010; % ~3.6 arcmin
	errd = 0.0010; % ~3.6 arcmin

	% Find list of sources to operate on
	cat1a = [cat1(ind1).alpha];
	cat1d = [cat1(ind1).delta];
	cat2a = [cat2.alpha];
	cat2d = [cat2.delta];

	 if (ep1 == true) % Cat1 is B1950
		if (ep2 == false) % Cat2 is J2000, the most likely scenario	
			% Convert the B1950 positions to J2000
			[x y z] = sph2cart (cat1a, cat1d, 1); 
			b1950cart = [x' y' z'];
			j2kcart = b1950cart * precMat;
			[cat1a2k, cat1d2k, ~] = ...
						cart2sph (j2kcart(:,1), j2kcart(:,2), j2kcart(:,3));
			cat1a2k (cat1a2k < 0) = cat1a2k (cat1a2k < 0) + 2*pi;
		else
			% Both catalogs are B1950
			cat1a2k = cat1a; cat1d2k = cat1d;
		end;
	else            % Cat1 is J2000
		if (ep2 == true)  % Cat2 is B1950
			[x y z] = sph2cart (cat2a, cat2d, 1); 
			b1950cart = [x' y' z'];
			j2kcart = b1950cart * precMat;
			[cat2a, cat2d, ~] = ...
						cart2sph (j2kcart(:,1), j2kcart(:,2), j2kcart(:,3));
			cat2a (cat2a < 0) = cat2a (cat2a < 0) + 2*pi;
			cat1a2k = cat1a; cat1d2k = cat1d; % For consistency in naming.
		else		% Both are J2000
			cat1a2k = cat1a; cat1d2k = cat1d;
		end;
	end;

	ind2 = zeros (1, length (ind1));
	% Carry out search for corresponding sources in both catalogs.
	for src = 1:length (ind1)
		% Make a list of sources with RA falling within chosen limits.
		rind = find ((cat2a < (cat1a2k(src)+errr)) & ...
					 (cat2a > (cat1a2k(src)-errr)));	

		if (length(rind) == 0)
			fprintf (1, '### No matching src found for source %d, based on RA. Error thresh: %f rad.\n', ...
					 cat1a2k(src), errr);
		else
			fprintf (1, '--> Src cat1 %d: %d candidates found, based on RA match.\n', ...
					 ind1(src), length(rind));
		end;

		% Check the dec of candidate sources
		for rsrc = 1:length (rind)		
			if ((cat2d (rind (rsrc))<(cat1d2k (src)+errd)) &  ...
				(cat2d (rind (rsrc))>(cat1d2k (src)-errd)))	
					ind2 (src) = rind(rsrc);
					if (debug > 0)
						fprintf (1, '  # Matched src cat2 %d. cat2d: %f, cat1d: %f\n', rind(rsrc), cat2d (rind(rsrc)), cat1d2k (src));
					end;
			end;
		end;
		if (ind2(src) == 0)
			fprintf (1, '### No matching src found for source %d, based on DEC. Error thresh: %f rad.\n', ...
					 cat1a2k(src), errd);
		end;
	end;
	
