% Script to generate statistics of 'causally' calibrated images, after 
% reading raw data from .bin files;
% NOTE: Expects the filenames to be defined in a cell variable 'fnames'.
% Pep/21May12

nfiles = length (fnames);
disp(' ');
disp (['Found ' num2str(nfiles) ' Files']);
figind = zeros (1, nfiles);
fid = zeros (1, nfiles);
acc = complex (zeros (288, 288, nfiles));
calvis = complex (zeros (288, 288, nfiles));
sigman = complex (zeros (288, 288, nfiles));
minarr = zeros (1, nfiles);
maxarr = minarr;
meanarr = minarr;
vararr = minarr;
% sigmas = zeros (1, 10, nfiles)); 
gainsol = complex (zeros (1, 288, nfiles));
saveasavi = 1;

tobs = zeros (1, nfiles);
fobs = zeros (1, nfiles);
duv = 2;
Nuv = 1000;                    % size of gridded visibility matrix
uvpad = 1024;                  % specifies if any padding needs to be added

% Local horizon based coordinates of array in ITRF
load ('poslocal.mat', 'posITRF', 'poslocal'); 

% Generate uv coordinates in local horizon coordinate system, needed for imaging
uloc = meshgrid (poslocal(:,1)) - meshgrid (poslocal (:,1)).';
vloc = meshgrid (poslocal(:,2)) - meshgrid (poslocal (:,2)).';
% normal to CS002 field (ITRF)
normal = [0.598753, 0.072099, 0.797682].'; 

for ind =1:nfiles
 fid(ind) = fopen (fnames{ind}, 'rb');
 outfname = [fnames{ind} '_stat.mat'];
 %if saveasavi == 1
 %  disp ('### NOTE: Saving movie without displaying!');
 %  disp ('### NOTE: Only saves a single file now!');
   % figind (ind) = figure ('Visible', 'Off');
   % aviobj  = avifile([fnames{ind} '.avi']);
 %else
 %  figind (ind) = figure();
 %end 
end
t=0;

% Calibrate from one timeslice of all files
for ind = 1:nfiles
  [acm(:,:,ind), tobs(:,ind), fobs(:,ind)] = readms2float (fid(ind), -1, -1);
  [thsrc_cat, phisrc_cat, thsrc_wsf, phisrc_wsf, calvis(:,:,ind), gainsol(:,:,ind), sigmas, sigman(:,:,ind), good(ind)] = pelican_sunAteamsub (acm(:,:,ind), tobs(ind), fobs(ind), eye(288), 2);
end 
   try
        while (feof (fid(1)) ~= 1 ) % Arbit. taking first file.
            t = t + 1;
            % Read in acm
	    for ind = 1:nfiles
  		[acm(:,:,ind), tobs(:,t), fobs(:,ind)] = readms2float (fid(ind), -1, -1);
                % Calibrate
                [thsrc_cat, phisrc_cat, thsrc_wsf, phisrc_wsf, calvis(:,:,ind), gainsol(:,:,ind), sigmas, sigman(:,:,ind), good(ind)] = pelican_sunAteamsub (acm(:,:,ind), tobs(t), fobs(ind), eye(288), 2);
                % [thcalvis, gainsol, sigmas, sigman, good] = pelican_sunAteamsub (acm , tobs, freq, eye(288), 2);
                calvis = (gainsol(:, ind)' * gainsol(:,ind)) .* acm(:,:,ind);
                % Generate per-subband images from each subband
                [calmap, ~] = fft_imager (calvis(:,:,ind), uloc, vloc, duv, Nuv, uvpad);
	 	minarr(t) = min(min(abs(calmap)));
	 	maxarr(t) = max(max(abs(calmap)));
		meanarr(t) = mean(mean(abs(calmap)));
		vararr(t) = var(var(abs(calmap)));
		disp ([num2str(tobs(:,t)) ' ' num2str(minarr(t)) ' ' num2str(maxarr(t)) ' ' num2str(meanarr(t)) ' ' num2str(vararr(t))]);
	   end
	end
    catch err
        disp ('Error encountered! Saving variables to disk..');
	save (outfname, 'tobs', 'minarr', 'maxarr');
        rethrow (err);
    end
    disp ('Saving to disk..');
    save (outfname, 'tobs', 'minarr', 'maxarr');
