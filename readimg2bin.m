% Script to read images written to binary files by wrimg2bin.m.
% pep/21Jan13

% Arguments:
%	fid: File id of the binary file to which calibrated visibilities should be
%		 written.

% Returns:
%	img.tobs: The MJD sec. time corresponding to this map.
%	img.freq: The Freq. in Hz. corresponding to this map.
% 	img.pix2laxis: Number of pixels along the l-axis.
% 	img.pix2maxis: Number of pixels along the m-axis.
%	img.map: Matrix containing the (square) image of size (pix2laxis X pix2maxis)
%			 to be written out. This is always real.
%   img.l,m : Vector containing the l and m coordinates corresponding to 
%				generated image.


function [img] = readimg2bin (fid, skip)

	if (fid < 0)
		disp ('readimg2bin: Inconsistent file id!');
		return;
	end;


	% Ideally, just an fseek should work, but recsize should be determined 
	% everytime. Doing the lazy thing for now...
	for ind = 0:skip  % Because skip == 0 should result in one record being read.
		img.tobs = fread (fid, 1, 'double');	
		if (isempty(img.tobs) == 1) 
			disp('readimg2bin: End of file reached!'); 
			return;
		end;

		img.freq      = fread (fid, 1, 'double');	
		img.pix2laxis = single (fread (fid, 1, 'float32'));
		img.l = single (fread (fid, img.pix2laxis, 'float32'));
		img.pix2maxis = single (fread (fid, 1, 'float32'));
		img.m = single (fread (fid, img.pix2maxis, 'float32'));
		map = single (fread (fid, img.pix2laxis * img.pix2maxis, 'float32'));
		img.map = reshape (map, img.pix2laxis, img.pix2maxis); 
	end;
