% Script to display images stored in .bin files, as generated by wrimg2bin.m
% pep/28Jan13
% Arguments:
%	fname  : image binary file name
%	colrng : caxis color range.
%	figarea: The area of the image to display, between 0 or 1.
%	offset : The record offset from which to start showing images.
%	skip   : Skip some number of images.
%	nrecs  : Number of images to show. -1 => show all images.
%	strail : Generate a star trail by accumulating images.
%	movie  : Flag to control creation of a .avi movie.
%			 -1 => Write out frames as png.
%			  0 => Just display on screen.
%			  1 => Write out frames as a single .avi
%
% Returns :
%          void

function showbinimages (fname, colrng, figarea, offset, skip, nrecs, strail, movie)
	fid = fopen (fname, 'rb');
	if (fid < 0)
		disp ('showbinimages: fid < 0! Quitting.');
		return;
	end;

	% Open movie related
	if (movie == 1)
		movname = strcat (fname, '.avi');
		vidObj = VideoWriter (movname, 'Indexed AVI');
		open (vidObj);
	end;

	if (movie == -1)
		mkdir ('binimg');
	end;

	img = readimg2bin (fid, 0);
	if (offset > 1)
		fprintf (2, 'Moving to image offset %d\n', offset);
		% NOTE that getting the size of a structured variable does not work like this...
		% tmp = whos ('img');
		% recsize = tmp.bytes;
		% NOTE: Hardcoded size of img constituents!
		recsize = 8 + ... % Frequency of image (double)
				  4 + ... % img.pix2laxis (float)
				  4*img.pix2laxis + ... % l axis values (floats)
				  4 + ... % img.pix2maxis (float)
				  4*img.pix2maxis + ... % m axis values (floats)
				  4*img.pix2maxis*img.pix2laxis; % Image pixel values (floats)

% THIS DID NOT WORK AND LED TO MOVING TO WRONG OFFSETS WITHIN THE FILE. pep/13May14
%		if (fseek (fid, (offset-1)*recsize, 'bof') < 0)
%			err = MException('showbinimages:fseek OutOfRange', ...
%        			'recoffset is outside expected range');	
%			error ('readms2float: seek error!');
%			throw (err);
%		end;
		for ind = 1:offset
			img = readimg2bin (fid, 0);
		end;
	end;
	map = reshape (img.map, img.pix2laxis, img.pix2maxis);
    mask = NaN(size (map));
	if (isempty(figarea))
		figarea = 1;
	end;
    mask(meshgrid(img.l).^2 + meshgrid(img.m).'.^2 < figarea) = 1;
	fprintf (1, sprintf ('Obs   info: Time: %.2f, freq: %.2f\n', ...
			 img.tobs, img.freq));
	fprintf (1,sprintf ('Image info: %fX%f, %f<l<%f, %f<m<%f\n',img.pix2laxis,...
			 img.pix2maxis, min(img.l), max(img.l), min(img.m), max(img.m)));


	if (nrecs < 0)
		fdir = dir (fname);
		filesize = fdir.bytes;
		imgwhos = whos ('img');
		imgsize = imgwhos.bytes;
		nrecs = filesize/imgsize;
	end;

	fig1 = figure;
	winsize = get (fig1, 'Position');
	winsize (1:2) = [0 0];

	for im = offset:offset+nrecs
		img = readimg2bin (fid, skip);
		if (strail == 1)
			map = (map + reshape (img.map, img.pix2laxis, img.pix2maxis))/2;
		else
			map = reshape (img.map, img.pix2laxis, img.pix2maxis);
		end;
    	fprintf(1, 'Showing image %03d of %03d, Time: %.2f\n', im, nrecs, ...
			  (img.tobs));
    	imagesc(img.l, img.m, real(map .* mask));
    
    	set(gca, 'FontSize', 16);
    	title (sprintf ('Time: %s, %.2f, Freq: %.2f', datestr(mjdsec2datenum(img.tobs)), img.tobs, img.freq));
	    axis equal
   		axis tight
   		set (gca, 'YDir', 'Normal'); % To match orientation with station images
   		set (gca, 'XDir', 'Reverse'); % To match orientation with station images

    	ylabel('South $\leftarrow$ m $\rightarrow$ North', 'interpreter', 'latex');
    	xlabel('East $\leftarrow$ l $\rightarrow$ West', 'interpreter', 'latex');
		if (~isempty (colrng))
			caxis (colrng);
		end;
    
    	set(colorbar, 'FontSize', 16);
		if (movie == 1)
			currFrame = getframe (fig1, winsize);
			writeVideo (vidObj, currFrame);
		end;

		if (movie == -1)
			currFrame = getframe (fig1, winsize);
			fimgname = sprintf ('binimg/%.0f.png', img.tobs);
			imwrite (currFrame.cdata, fimgname, 'png');
		end;
		pause (0.01);
		% pause;
	end;

	if (movie == 1)
		close (vidObj);
	end;
	fclose (fid);

	 
