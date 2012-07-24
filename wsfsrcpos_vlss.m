% Program to estimate the positions of dominant sources in a visibility set
% using Weighted Subspace Fitting method. The initial model is now
% generated using the positions of the dominant VLSS sources.
% pep/22Mar12

function wsfsrcpos_vlss (fname)
    global stop_processing
    %% Basic initialization
    freq = 59951782.2;
    C = 299792458; % speed of light
    l = [-1:0.005:1]; m = l;
    lon = 6.869837540;                         % longitude of CS002 in degrees
    lat = 52.915122495;                        % latitude of CS002 in degrees
    restriction = 4;             % Avoid uv values with <10 wavelengths length
    maxrestriction = 60;          % Avoid uv values with >60 wavelengths length
    % fnames = {'SB000_1415-1420.bin', '1222_1227.bin'};
    nfiles = 2;
    Nelem = 288; 
    nblines = Nelem * (Nelem + 1)/2; 
    load srclist3CR
    load ('poslocal.mat', 'posITRF', 'poslocal'); 
    normal = [0.598753, 0.072099, 0.797682].';
    srcsel =  [324, 283, 88, 179, 0]; % A team from 3CR catalog
    % srcsel = [57759, 66461, 66462, 16848, 39974, 0]; % Ateam from VLSS catalog
    maxtslots = 200;

    % VLSS reading in section
    nvlsssrc = 10; % Take the first 10 sources with highest flux
    vlss = fitsread ('CATALOG.FIT', 'BinTable');

    % [tmp ind] = sort (vlss{3}, 'descend');
    % srcvlsssel = ind(1:nvlsssrc);
    % disp (['Total flux in top' num2str(nvlsssrc) ' sources (Jy): '
    % num2str(sum(vlss{3}(srcvlsssel)))]);
    for ind = 1:length (srcsel)-1
       vlsscat(ind).alpha = vlss{1}(srcsel(ind))*pi/180; % convert to rad.
       vlsscat(ind).delta = vlss{2}(srcsel(ind))*pi/180;
       vlsscat(ind).flux  = vlss{3}(srcsel(ind));
    end
    vlsscat(length(srcsel)).alpha = 0;
    vlsscat(length(srcsel)).delta = 0;
    vlsscat(length(srcsel)).flux = 0;
    % NOTE: Setting srcsel to consecutive indices due to collation in vlsscat
    % srcsel = [1:length(srcsel)-1 0];
    catalog = srclist3CR;


    fid = fopen (fname, 'rb');
    outfilename = strrep (fname, '.bin', '_wsflmpos.mat');
    disp(['Writing estimated  and catalog positions to file: ' outfilename]);

    t_obs = fread (fid, 1, 'double');
    disp (['Timeslice: ' num2str(t_obs, '%f') ' (MJD secs) for initial estimation']);
    

    % Reading real and imaginary, available as a stream of floats.
    % even floats being real parts, odd floats being imag
    a = fread (fid, 2*nblines, 'float'); % Read one ccm worth
    comp = complex (a(1:2:length(a)), a(2:2:length(a))); % to complex
    % create instantaneous ccm from vector
    acm = triu (ones (Nelem));
    acm (acm == 1) = comp;
    acc = acm + acm' - diag(diag(acm));
    
    t_obs1 = fread (fid, 1, 'double');
    tslices2cal = floor ((2*60)/(t_obs1 - t_obs)); % recalibrate every 2 mins 
    disp (['Time resolution of observations: ' num2str(t_obs1 - t_obs) ' secs']);
    disp (['Calibrating every ' num2str(tslices2cal) ' timeslices']);
    
    % Convert back to Julian day, as desired by track_statcal. NOTE:
    % JulianDay () should not be called now!
    t = 1;
    t_obs_store (t) = t_obs;
    t_obs  = t_obs/86400 + 2400000.5; % Convert to day units
    t_obs1 = t_obs1/86400 + 2400000.5;
   
    % Go to the beginning of the file
    fseek (fid, 0, -1);

    %% calibrate using source positions from catalog and a baseline
    % restriction of 10 wavelengths or 60 meters (whichever is smaller)
    % This results in estimated source fluxes of the model sources
    [cal1, sigmas1, Sigman1] = statcal_vlss(acc, t_obs, freq, posITRF, srcsel, normal, 10, 60, eye(Nelem), catalog);
    cal_store (:,1) = cal1;
    sigmas_store = {};
    sigmas_store {1} = sigmas1;
    sigman_store(:,:,1) = Sigman1;
    
    sel = sigmas1 > 0.001;
    nsrc = sum(sel);
    calrun = 1;
    stop_processing = 0 ;
    %%
    % stop_proc ();
    try
        while (feof (fid) ~= 1 & stop_processing == 0)
            %keypress = isfigure (f);
        % for i=1:10

            % disp (['  --> Estimating source positions from data...']);
        %     sel = sigmas1 > 0.01;
        %     nsrc = sum(sel);
            % use source positions from 3C catalog as initial estimate
            % load srclist3CR
            if (sum(srcsel == 0) ~= 0)
                [raSun, decSun] = SunRaDec(t_obs);                % J2000
                rasrc = [catalog(srcsel(srcsel ~= 0)).alpha, raSun].';   % B1950
                decsrc = [catalog(srcsel(srcsel ~= 0)).delta, decSun].'; % B1950
            else
                rasrc = [catalog(srcsel(srcsel ~= 0)).alpha];  % B1950
                decsrc = [catalog(srcsel(srcsel ~= 0)).delta]; % B1950
            end
            [lsrc0, msrc0] = radectolm(rasrc(sel), decsrc(sel), t_obs, lon, lat);
            lsrc_cat (:,t) = lsrc0;
            msrc_cat (:,t) = msrc0;
            


            %% Apply calibration on data
            % acccal = (cal1'*cal1).*acc; 
            acccal = acc;
            %  find source positions using WSF algorithm
            [v, d] = eig(acccal);
            [d, order] = sort(diag(abs(d)), 'descend');
            v = v(:, order);
            Wopt = (diag(d(1:nsrc)) - mean(diag(squeeze(Sigman1))) * eye(nsrc))^2 / diag(d(1:nsrc));
            Es = v(:, 1:nsrc);
            EsWEs = Es * Wopt * Es';
            theta = fminsearch(@(x) WSFcostfun(x, EsWEs, diag(conj(1 ./ cal1)), freq, poslocal(:, 1).', poslocal(:, 2).'), [lsrc0; msrc0]);
            lsrchat = theta(1:nsrc);
            msrchat = theta(nsrc+1:end);
            lsrc_wsf (:,t) = lsrchat;
            msrc_wsf (:,t) = msrchat;
            t = t + 1;

            t_obs = fread (fid, 1, 'double');
            if (isempty (t_obs))
                break;
            end;
            disp (['Timeslice: ' num2str(t) ' (' num2str(t_obs, '%f') ')']);
            % Convert back to Julian day, as desired by track_statcal. NOTE:
            % JulianDay () should not be called now!
            t_obs_store (t) = t_obs;
            t_obs = t_obs/86400 + 2400000.5;

            % Reading real and imaginary, available as a stream of floats.
            % even floats being real parts, odd floats being imag
            a = fread (fid, 2*nblines, 'float'); % Read one ccm worth
            if (isempty (a))
                break;
            end;
            comp = complex (a(1:2:length(a)), a(2:2:length(a))); % to complex
            % create instantaneous ccm from vector
            acm = triu (ones (Nelem));
            acm (acm == 1) = comp;
            acc = acm + acm' - diag(diag(acm));
% 
%             for t_ind = 1:50
%                 t_obs = fread (fid, 1, 'double');
%                 if (isempty(t_obs))
%                     break;
%                 end
%                 a = fread (fid, 2*nblines, 'float'); % Read one ccm worth
%                 if (isempty(a))
%                     break;
%                 end
%             end    
            if mod (t, tslices2cal) == 0
                [cal1, sigmas1, Sigman1] = statcal_vlss(acc, t_obs, freq, posITRF, srcsel, normal, 10, 60, eye(Nelem), catalog);
                calrun = calrun + 1;
                cal_store (:,calrun) = cal1;
                sigmas_store {calrun} = sigmas1;
                sigman_store (:,:,calrun) = Sigman1;

            end    
        end
    catch err
        disp ('There was an error! Saving current set of variables...');
        save (outfilename, 'lsrc_cat', 'msrc_cat', 'lsrc_wsf', 'msrc_wsf', 'cal_store', 'sigmas_store', 'sigman_store', 't_obs_store');
        rethrow (err);
    end   
    disp ('WSF Estimation completed. Storing variables to file');
    save (outfilename, 'lsrc_cat', 'msrc_cat', 'lsrc_wsf', 'msrc_wsf', 'cal_store', 'sigmas_store', 'sigman_store', 't_obs_store');