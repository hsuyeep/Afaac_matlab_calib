% Script to write out generated images as floats matrices, in the format accepted
% by readimg2bin.m
% pep/21Jan13

% Arguments:
%	fid: File id of the binary file to which calibrated visibilities should be
%		 written.
%	img.tobs: The MJD sec. time corresponding to this map.
%	img.freq: The Freq. in Hz. corresponding to this map.
% 	img.pix2laxis: Number of pixels along the l-axis.
% 	img.pix2maxis: Number of pixels along the m-axis.
%	img.map: Matrix containing the (square) image of size (pix2laxis X pix2maxis)
%			 to be written out. This is always real.
%   img.l,m : Vector containing the l and m coordinates corresponding to 
%				generated image.


function wrimg2bin (fid, img)

	if (fid < 0)
		disp ('wrimg2bin: Inconsistent file id!');
		return;
	end;

	if ((img.pix2laxis ~= length (img.l)) || (img.pix2maxis ~= length (img.m)))
		disp ('wrimg2bin: Inconsistent number of pixels!');
		return;
	end;

	im = whos ('img');
    recsize = im.bytes;

	disp (sprintf('wrimg2bin: Writing tobs: %f, Freq: %f', img.tobs, img.freq));
	count = fwrite (fid, img.tobs, 'double');
	count = fwrite (fid, img.freq, 'double');
	count = fwrite (fid, single(img.pix2laxis), 'float32');
	count = fwrite (fid, single(img.l), 'float32');
	count = fwrite (fid, single(img.pix2maxis), 'float32');
	count = fwrite (fid, single(img.m), 'float32');
	count = fwrite (fid, single(img.map(:)), 'float32');
