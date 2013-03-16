% Script to generate closure quantities of calibrated visibilities to check
% accuracy of calibration.
% NOTE: Closure phases are expected to be close to 0 only for a single point 
% source in the field! Not for a general sky...
% pep/03Dec12

function genclos(fname, nrecs)
 	if nrecs == -1
		[nrecs, tmin, tmax, dt] = getnrecs (fname);
	end;
	col = {'bo', 'mo', 'ro', 'ko', 'go', 'yo', 'wo', 'co'};
	fid = fopen (fname, 'rb');
	[acc, tobs, freq] = readms2float (fid, -1, -1, 288);
	acc (isnan(acc)) = 0; % NOTE NOTE! Norm fails if members are NaNs!
	Nelem = size (acc, 1);
	cnt = 1;
	% ant_set = [23 + 1:43:288]; % Smaller set of antennas to work with.
	ant_set = [23 + 1:80:288]; % Smaller set of antennas to work with.
	nants = length (ant_set);
	ntrips = nants * (nants - 1) * (nants - 2) / 6;
	fprintf (1, 'genclos: %d antennas, %d triples for each of %d recs.\n',...
			 nants, ntrips, nrecs);
	clos_amp = zeros (ntrips, nrecs);
	clos_ph = clos_amp;
	

	% NOTE: Taking visibility triples from all baselines results in 
	% n(n-1)(2n-1)/6 pairs ~ 8million combinations, so taking representative
	% baselines from eahc station with every other station. 
	for ts=1:nrecs
%	  	for ind = 1:Nelem % Set of all antennas
%     		for jind = ind+1:Nelem
%         		for kind = jind+1:Nelem

		cnt = 0;
		for ind = 1:nants-2
			for jind = ind+1:nants-1
          		for kind = jind+1:nants
			
%             		if (flag (ind, jind) == 0 || flag (jind, kind) == 0 ... 
%	 						|| flag (kind, ind) == 0)
%                 		flag (ind ,jind) = 1; flag (jind, ind) = 1;
%                 		flag (jind,kind) = 1; flag (kind,jind) = 1;
%                 		flag (kind, ind) = 1; flag (ind ,kind) = 1;
                 	cnt = cnt + 1;
                 	clos_amp (cnt, ts) =  abs  (acc(ind,jind)) ... 
	 						* abs  (acc(jind,kind)) * abs  (acc(kind,ind));
	                clos_ph  (cnt, ts) = angle(acc(ind,jind)) ...
	 					+ angle(acc(jind,kind)) + angle(acc(kind,ind));
	 				% clos_ph (cnt) = sin (clos_ph(cnt));

             	end
         	end
     	 end
		[acc, tobs, freq] = readms2float (fid, -1, -1, 288);
		acc (isnan(acc)) = 0; % NOTE NOTE! Norm fails if members are NaNs!
	end

	% disp (['Total triples: ' num2str(cnt)]);
	subplot (3,1,1);
	for ind = 1:min (length(col), ntrips)
		plot (clos_ph(ind,:), char(col(ind)));
		hold on;
	end;
	if (ntrips > length(col))
		for ind = length(col):ntrips
			plot (clos_ph(ind,:), char(col(end)));
			hold on;
		end;
	end;

	title (sprintf ('Closure phase for %d triples from %d antennas.\n', ...
			ntrips,  nants));

	subplot (3,1,2);
	for ind = 1:min (length(col), ntrips)
		plot (clos_amp(ind,:), char(col(ind)));
		hold on;
	end;
	if (ntrips > length(col))
		for ind = length(col):ntrips
			plot (clos_amp(ind,:), char(col(end)));
			hold on;
		end;
	end;
	title (sprintf ('Closure amplitude for %d triples from %d antennas.\n', ...
			ntrips, nants));

	subplot (3,1,3);
	imagesc (hist (clos_ph', 100));

	% normts(1) = norm (acc, 'fro');
	% for ts=2:nrecs
	% end;
