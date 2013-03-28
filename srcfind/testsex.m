% Script to test the source extraction using gaussian fitting.
% pep/02Jan13
clear all; close all;
addpath '../'
debug = 2;
radec = 0;
duv = 2.5;						% Default, reassigned from freq. of obs. to
								% image just the full Fov (-1<l<1)
Nuv = 500; %1000                % size of gridded visibility matrix
uvpad = 512; %1024              % specifies if any padding needs to be added
nfacet = 3; facetsize = 256;

% Local horizon based coordinates of array in ITRF
load ('../poslocal.mat', 'posITRF', 'poslocal'); 
uloc = meshgrid (poslocal(:,1)) - meshgrid (poslocal (:,1)).';
vloc = meshgrid (poslocal(:,2)) - meshgrid (poslocal (:,2)).';

% location of calibrated visibilities:  Daytime data
% fname = '~/WORK/AARTFAAC/rawdat/SB000_ch30-35_5sec_3hr_cal.bin'; 

% location of calibrated visibilities:  Nighttime data
fname = '~/WORK/AARTFAAC/Reobs/11Jul12/LBA_OUTER_BAND60/SB001_LBA_OUTER_4min_1ch_cal.bin';
fid = fopen (fname, 'rb');

% longitude and latitude of CS002 (AARTFAAC center).
lon = 6.869837540;                         % longitude of CS002 in degrees
lat = 52.915122495;                        % latitude of CS002 in degrees 


[acc, tobs, freq] = readms2float (fid, -1, -1, 288);
tobs_jd = tobs/86400. + 2400000.5; % Convert to JD.
lambda = 299792458/freq; 		% in m.
duv = lambda/2;
dl = (299792458/(freq * uvpad * duv)); % dimensionless, in dir. cos. units
lmax = dl * uvpad / 2;

% Check if Sun subtraction needs to be carried out.
% Locate the Sun in this image.
[raSun, decSun] = SunRaDec(tobs_jd);
[l_sun, m_sun] = radectolm(raSun, decSun, tobs_jd, lon, lat);
% TODO: Check if the sun is above the horizon at this time.	
% TODO: Generate Sun subtracted visibilities.

% Generate a map from calibrated visibilities
[calmap, l, m] = genmosaic(acc, tobs, freq, nfacet, facetsize, uloc(:), vloc(:), posITRF, 0);

% Get image statistics
[dr, sig] = getimagedr (calmap, 32, 3); % 64-pixel square, 3sigma threshold.

extract (calmap, l, m, tobs, freq, 300, sig, 5, debug);

% imagesc (l, m, calmap);
