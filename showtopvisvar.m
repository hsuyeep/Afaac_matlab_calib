% Script to show the timeseries of visibilities and their variances,
% As extracted from the TOP MS file.
% pep/17Dec12

fid = fopen ('~/WORK/TOP_CygA/L77781_SAP000_SB123_uv.wbin', 'rb');
nchan = 3; nant = 34; % Hardcoded

for ind = 1:100
	for ch=1:3
		[acc, weight, tobs, freq] = readms2weightfl(fid, -1, -1, nant);
		subplot (3,2,1);
		imagesc (20*log10(abs(acc)));
		colorbar;
		title (sprintf ('Visibility. Time: %f, freq: %f', tobs, freq));
		
		subplot (3,2,2);
		imagesc (weight);
		colorbar;
		title (sprintf ('Weights. Time: %f, freq: %f', tobs, freq));

		[acc, weight, tobs, freq] = readms2weightfl(fid, -1, -1, nant);
		subplot (3,2,3);
		imagesc (20*log10(abs(acc)));
		colorbar;
		title (sprintf ('Visibility. Time: %f, freq: %f', tobs, freq));
		
		subplot (3,2,4);
		imagesc (weight);
		colorbar;
		title (sprintf ('Weights. Time: %f, freq: %f', tobs, freq));

		[acc, weight, tobs, freq] = readms2weightfl(fid, -1, -1, nant);
		subplot (3,2,5);
		imagesc (20*log10(abs(acc)));
		colorbar;
		title (sprintf ('Visibility. Time: %f, freq: %f', tobs, freq));
		
		subplot (3,2,6);
		imagesc (weight);
		colorbar;
		title (sprintf ('Weights. Time: %f, freq: %f', tobs, freq));
		% pause;
	end;
end;
