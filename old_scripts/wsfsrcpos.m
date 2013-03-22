% Program to estimate the positions of dominant sources in a visibility set
% using SJW's Weighted Subspace Fitting method, now operating using ITRF.
% pep/16Mar12

%% Basic initialization
freq = 59951782.2;
C = 299792458; % speed of light
l = [-1:0.005:1]; m = l;
lon = 6.869837540;                         % longitude of CS002 in degrees
lat = 52.915122495;                        % latitude of CS002 in degrees
restriction = 4;             % Avoid uv values with <10 wavelengths length
maxrestriction = 60;          % Avoid uv values with >60 wavelengths length
fname = 'SB000_ch30-35_5sec_3hr.bin';
% fname = 'SB000_1415-1420.bin';
skiprecs = 120;
nfiles = 2;
Nelem = 288; 
nblines = Nelem * (Nelem + 1)/2; 
load srclist3CR
load ('poslocal.mat', 'posITRF', 'poslocal'); 
normal = [0.598753, 0.072099, 0.797682].';
srcsel =  [324, 283, 88, 179, 0];

fid = fopen (fname, 'rb');
outfilename = strrep (fname, '.bin', '_wsfthphipos.mat');
disp(['Writing estimated  and catalog positions to file: ' outfilename]);

t_obs1 = fread (fid, 1, 'double');
disp (['Timeslice: ' num2str(t_obs1, '%f') ' (MJD secs) for initial estimation']);
% Convert back to Julian day, as desired by track_statcal. NOTE:
% JulianDay () should not be called now!
t_obs = t_obs1/86400 + 2400000.5;

% Reading real and imaginary, available as a stream of floats.
% even floats being real parts, odd floats being imag
a = fread (fid, 2*nblines, 'float'); % Read one ccm worth
comp = complex (a(1:2:length(a)), a(2:2:length(a))); % to complex
% create instantaneous ccm from vector
acm = triu (ones (Nelem));
acm (acm == 1) = comp;
acc = acm + acm' - diag(diag(acm));


%% calibrate using source positions from 3C catalog and a baseline
% restriction of 10 wavelengths or 60 meters (whichever is smaller)
% This results in estimated source fluxes of the model sources
% [cal1, sigmas1, Sigman1] = statcal(acc, t_obs, freq, posITRF, srcsel, normal, 10, 60, eye(Nelem));
[cal1, sigmas1, Sigman1] = statcal_stefcal (acc, t_obs, freq, posITRF, srcsel, normal, 10, 60, eye(Nelem));

% whitening of the array covariance matrix (required for DOA estimation)
% NOTE: Whitening before calibration causes calibration to fail!
acc = acc ./ sqrt(diag(acc) * diag(acc).');
acc(isnan(acc)) = 0; % Get rid of NaNs and infinities which kill eig ()
acc(isinf(acc)) = 0;

disp('First round of antenna calibration completed, sigmas1:'); disp (sigmas1');
disp(['NOTE: Skipping ' num2str(skiprecs) 'Timeslices!']);
sel = sigmas1 > 0.01;
nsrc = sum(sel);
disp (['WSF: Working with ' num2str(nsrc) ' selected sources, t_obs: ' num2str(t_obs)]);
% lsrc_cat = zeros (nsrc, 200);
% msrc_cat = lsrc_cat;
% lsrc_wsf = lsrc_cat;
% msrc_wsf = lsrc_cat;
first = 1;
t = 1;
%%
while (feof (fid) ~= 1)
% for i=1:10
    
    % disp (['  --> Estimating source positions from data...']);
%     sel = sigmas1 > 0.01;
%     nsrc = sum(sel);
    % use source positions from 3C catalog as initial estimate
    % load srclist3CR
    if (sum(srcsel == 0) ~= 0)
        [raSun, decSun] = SunRaDec(t_obs);                % J2000
        rasrc = [srclist3CR(srcsel(srcsel ~= 0)).alpha, raSun].';   % B1950
        decsrc = [srclist3CR(srcsel(srcsel ~= 0)).delta, decSun].'; % B1950
        epoch = [true(length(srcsel(srcsel ~= 0)), 1); false];
    else
        rasrc = [srclist3CR(srcsel(srcsel ~= 0)).alpha];  % B1950
        decsrc = [srclist3CR(srcsel(srcsel ~= 0)).delta]; % B1950
        epoch = true(length(srcsel), 1);
    end
    srcpos0 = radectoITRF(rasrc(sel), decsrc(sel), epoch(sel), t_obs);
    thsrc0 = asin(srcpos0(:, 3));
    phisrc0 = atan2(srcpos0(:, 2), srcpos0(:, 1));
    % [lsrc0, msrc0] = radectolm(rasrc(sel), decsrc(sel), t_obs, lon, lat);
    % lsrc_cat (:,t) = lsrc0;
    % msrc_cat (:,t) = msrc0;
    % lsrc_cat (:,t) = thsrc0;
    % msrc_cat (:,t) = phisrc0;
    
    
    %%  find source positions using WSF algorithm
    [v, d] = eig(acc);
    [d, order] = sort(diag(abs(d)), 'descend');
    v = v(:, order);
    Wopt = (diag(d(1:nsrc)) - mean(diag(squeeze(Sigman1))) * eye(nsrc))^2 / diag(d(1:nsrc));
    Es = v(:, 1:nsrc);
    EsWEs = Es * Wopt * Es';
    % theta = fminsearch(@(x) WSFcostfun(x, EsWEs, diag(conj(1 ./ cal1)), freq, poslocal(:, 1).', poslocal(:, 2).'), [lsrc0; msrc0]);
    % lsrchat = theta(1:nsrc);
    % msrchat = theta(nsrc+1:end);
    % lsrc_wsf (:,t) = lsrchat;
    % msrc_wsf (:,t) = msrchat;
    theta = fminsearch(@(x) WSFcostfunITRF(x, EsWEs, diag(conj(1 ./ cal1)), freq, posITRF), [phisrc0; thsrc0]);
    phisrchat = theta(1:nsrc);
    thsrchat = theta(nsrc+1:end);
    srcposhat = [cos(phisrchat) .* cos(thsrchat), sin(phisrchat) .* cos(thsrchat), sin(thsrchat)];
    disp('Position estimation using WSF done: ');
    disp ('Catalog positions: '); disp ([thsrc0 phisrc0]');
    disp ('WSF     positions: '); disp ([thsrchat phisrchat]');
    if first == 1
     first = 0;
     thsrc_cat  (:,1) = thsrc0;
     phisrc_cat (:,1) = phisrc0;
     thsrc_wsf  (:,1) = thsrchat;
     phisrc_wsf (:,1) = phisrchat;

     % srcpos_cat = {};
     % srcpos_cat{1} = srcpos0;
     % srcpos_wsf = {};
     % srcpos_wsf{1} = srcposhat;
     trec (1) = t_obs1;
    else 
     thsrc_cat  (:,t) = thsrc0;
     phisrc_cat (:,t) = phisrc0;
     thsrc_wsf  (:,t) = thsrchat;
     phisrc_wsf (:,t) = phisrchat;

     % srcpos_cat {t} = srcpos0;
     % srcpos_wsf {t} = srcposhat;
      trec(t) = t_obs1;
    end
    t = t + 1;
    
    % Read data for the next timeslice
    for i=1:skiprecs
       t_obs1 = fread (fid, 1, 'double');
       if isempty(t_obs1)
         break;
       end;
       a = fread (fid, 2*nblines, 'float'); % Read one ccm worth
       if isempty(a)
         break;
       end;
    end
    disp (['Timeslice: ' num2str(t) ' (' num2str(t_obs1, '%f') ')']);
    % Convert back to Julian day, as desired by track_statcal. NOTE:
    % JulianDay () should not be called now!
    t_obs = t_obs1/86400 + 2400000.5;

    % Reading real and imaginary, available as a stream of floats.
    % even floats being real parts, odd floats being imag
    comp = complex (a(1:2:length(a)), a(2:2:length(a))); % to complex
    % create instantaneous ccm from vector
    acm = triu (ones (Nelem));
    acm (acm == 1) = comp;
    acc = acm + acm' - diag(diag(acm));
    % whitening of the array covariance matrix (required for DOA estimation)
    acc = acc ./ sqrt(diag(acc) * diag(acc).');
    acc(isnan(acc)) = 0; % Get rid of NaNs
    acc(isinf(acc)) = 0;
end

% save (outfilename, 'lsrc_cat', 'msrc_cat', 'lsrc_wsf', 'msrc_wsf');
save (outfilename, 'thsrc_cat', 'phisrc_cat', 'thsrc_wsf', 'phisrc_wsf', 'trec');
