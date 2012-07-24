% Script to generate a series of images from a ms2float binary file
% NOTE: Expects the filenames to be defined in a cell variable 'fnames'.
% Pep/19May12
% Added ability to see multiple files (assumed to have different frequencies)
% in separate figure windows, as well as to save the images as an AVI file
% without displaying them.
% Pep/21May12

nfiles = length (fnames);
disp(' ');
disp (['Found ' num2str(nfiles) ' Files']);
figind = zeros (1, nfiles);
fid = zeros (1, nfiles);
acc = complex (zeros (288, 288, nfiles));
calvis = complex (zeros (288, 288, nfiles));
sigman = complex (zeros (288, 288, nfiles));
% sigmas = zeros (1, 10, nfiles)); 
gainsol = complex (zeros (1, 288, nfiles));
saveasavi = 0;

tobs = zeros (1, nfiles);
fobs = zeros (1, nfiles);
duv = 2;
Nuv = 500;                    % size of gridded visibility matrix
uvpad = 512;                  % specifies if any padding needs to be added
% Nuv = 1000;                    % size of gridded visibility matrix
% uvpad = 1024;                  % specifies if any padding needs to be added

% Local horizon based coordinates of array in ITRF
load ('poslocal.mat', 'posITRF', 'poslocal'); 

% Generate uv coordinates in local horizon coordinate system, needed for imaging
uloc = meshgrid (poslocal(:,1)) - meshgrid (poslocal (:,1)).';
vloc = meshgrid (poslocal(:,2)) - meshgrid (poslocal (:,2)).';
% normal to CS002 field (ITRF)
normal = [0.598753, 0.072099, 0.797682].'; 

for ind =1:nfiles
 fid(ind) = fopen (fnames{ind}, 'rb');
 if saveasavi == 1
   disp ('### NOTE: Saving movie without displaying!');
   disp ('### NOTE: Only saves a single file now!');
   figind (ind) = figure ('Visible', 'Off');
   aviobj  = avifile([fnames{ind} '.avi']);
 else
   figind (ind) = figure();
 end 
end
t=1;

% Calibrate from one timeslice of all files
for ind = 1:nfiles
  [acm(:,:,ind), tobs(:,ind), fobs(:,ind)] = readms2float (fid(ind), t, -1);
  [thsrc_cat, phisrc_cat, thsrc_wsf, phisrc_wsf, calvis(:,:,ind), gainsol(:,:,ind), sigmas, sigman(:,:,ind), good(ind)] = pelican_sunAteamsub (acm(:,:,ind), tobs(ind), fobs(ind), eye(288), 2);
end 
   try
        while (feof (fid(1)) ~= 1 ) % Arbit. taking first file.
            t = t + 1;
            % Read in acm
	    for ind = 1:nfiles
  		[acm(:,:,ind), tobs(:,ind), fobs(:,ind)] = readms2float (fid(ind), -1, -1);
                % Calibrate
                %[calvis, gainsol, sigmas, sigman, good] = pelican_sunAteamsub (acm , tobs, freq, eye(288), 2);
                calvis = (gainsol(:, ind)' * gainsol(:,ind)) .* acm(:,:,ind);
                % Generate per-subband images from each subband
                [calmap, ~] = fft_imager (calvis, uloc, vloc, duv, Nuv, uvpad);
		if saveasavi ~= 1
		  figure (figind(ind)); % Make current window visible, on top
                else
	          % Make the figure h current, but do not change its visibility
	          set (0, 'CurrentFigure', figind(ind));
                end
	        lab  = sprintf ('Time: %15.3f, Freq: %d',tobs(ind), int64(fobs(ind)));
	        imagesc (abs(calmap));
		colorbar;
	        title (lab);
		if saveasavi == 1
		  F = getframe(figind(ind));
                  aviobj = addframe(aviobj,F);
		end;
	   end
	end
    catch err
        disp ('Error encountered! Saving variables to disk..');
	if saveasavi == 1
          aviobj = close (aviobj);
        end 
        % for i=1:nfiles
        %   aviobj(ind) = close (aviobj(ind));
        % end
        rethrow (err);
    end
    if saveasavi == 1
      aviobj = close (aviobj);
    end 
    % for i=1:nfiles
    %    aviobj(ind) = close (aviobj(ind));
    % end
