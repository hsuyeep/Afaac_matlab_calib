% script to readin the VLSS catalog into a structure similar to the one
% used by SJW for representing the 3CR catalog.

% Arguments: 
%	fname  : Filename of the VLSS ascii catalog
% srcfluxlimit: All sources above the specified srcfluxlimit will be extracted.
% Returns :
%   srclistVLSS: A structure containing the extracted sources with names, 
%				 positions and fluxes.

% pep/20Jan12


function srclistVLSS = readVLSS (fname, srcfluxlimit)
	% srcnumlimit = 100;% Only the first 'srcnumlimit' sources are picked

	% In Jy. only sources having a flux above this are picked.
	if srcfluxlimit == 0
		srcfluxlimit = 20; 
	end;
	fvlss = fopen (fname, 'r'); 
	
	% Eliminate text in header, but print out relevant portions
	str = fgets (fvlss);
	while strcmp (strtok(str), 'RA(2000)') == 0
	  str = fgets (fvlss);
	  fprintf (1, '%s', str); 
	end
	
	% eliminate  text line containing units
	str = fgets (fvlss);
	
	srccnt = 1;
	% Not imposing a srcnumlimit
	% while (feof (fvlss) == 0 && srccnt < srcnumlimit)
	while (feof (fvlss) == 0)
		str = fgets (fvlss);
		[val cnt] = sscanf (str, '%d %d %f %d %d %f %f', 7);
		
		if val(7) > srcfluxlimit
			srclistVLSS (srccnt).name = 'test';
			% Convert positions to decimal hours, and then to radians.
			srclistVLSS (srccnt).alpha = (val(1) + val(2)/60 + val(3)/3600) ...
										 /(2*pi);
			if val (4) < 0
				srclistVLSS(srccnt).delta=(val(4)-val(5)/60-val(6)/3600)/(2*pi);
			else
				srclistVLSS(srccnt).delta=(val(4)+val(5)/60+val(6)/3600)/(2*pi);
			end 
			srclistVLSS(srccnt).flux = val (7);
			srccnt = srccnt + 1;
		end 
		
		str = fgets (fvlss);    % ignore next line which only defines errors
	end

	foutname =  sprintf ('srclistVLSS_%d.mat', srcfluxlimit);
	fprintf (1, 'Sources extracted: %d\n...Now writing to file %s.\n', ...
			 srccnt,foutname);
	save (foutname, 'srclistVLSS');
