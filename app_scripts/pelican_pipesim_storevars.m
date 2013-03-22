% Program to simulate a pelican pipeline, using the pelican_calib ()
% and fft_imager_sjw () functions.
% pep/23Mar12
% Updated to use readms2float.m

function pelican_pipesim_storevars (fname)
   %%
      Nelem = 288; 
    nblines = Nelem * (Nelem + 1)/2; 
    uvflag = eye (288); % Flag only the autocorrelations
    mins2cal = 2;
    debuglev = 2;
    skiprecs = 10;   % NOTE: Should not be 0, else same timeslice is processed!
    recoffset = 10;

    % FFT Imaging related
    % duv = 600/511;                % In meters, grid spacing in fft imager
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

    
    fid = fopen (fname, 'rb');
    [acm, tobs, fobs] = readms2float (fid, recoffset, -1);

   t = 0;
   try
        while (feof (fid) ~= 1 )
        % for i=1:5

            t = t + 1;
            if (isempty (tobs))
                break;
            end
            disp ('---->');
            outfilename = sprintf ('%10.0f_cal.mat',tobs);
            disp (['Timeslice: ' num2str(t) ' (' num2str(tobs, '%f') '), to file: ' outfilename]);

            
            % Calibrate this timeslice
            % [calvis, gainsol, sigmas, sigman, good] = pelican_calib (acc, tobs, freq, uvflag);
            [thsrc_cat, phisrc_cat, thsrc_wsf, phisrc_wsf, calvis, gainsol, sigmas, sigman, good] = pelican_sunAteamsub (acm, tobs, fobs, uvflag, 2);
            % Image this data
            % [calmap, ~] = fft_imager (calvis, uloc, vloc, duv, Nuv, uvpad);
            

            save (outfilename, 'calvis', 'gainsol', 'sigmas', 'sigman', 'tobs', 'fobs', 'good', 'thsrc_cat', 'phisrc_cat', 'thsrc_wsf', 'phisrc_wsf');
            [acm, tobs, fobs] = readms2float (fid, recoffset+t*skiprecs, -1);
        end
    catch err
        disp ('Error encountered! Saving variables to disk..');
        % save (outfilename, 'gainsol', 'sigmas', 'sigman', 'tobs', 'calmap', 'good');
        rethrow (err);
    end

