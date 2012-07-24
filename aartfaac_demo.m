% Demonstration of AARTFAAC imaging
% Carries out a basic calibration (via the Station calibration routines)
% of the superterp ACM timeseries made available as  matlab matrix files.
%  Uses the DFT imager to generate images, which are saved for every time
% slice.
% SJW, 11 January 2012

%% start with a clean workspace
clear all
close all

%% Extract out the timestamp of the observation and load cross correlation
%% matrix
ccm_fname = '../cookdat/sterp_tseries/SB000_133908_sterp_clean.mat';
load (ccm_fname);

[tok ccm_fname] = strtok (ccm_fname, '/');
a = size (ccm_fname);
while a(2) ~= 0
    [tok ccm_fname] = strtok (ccm_fname, '/');
    a = size (ccm_fname);
end
    
f_tobs = sscanf (tok, 'SB000_%2d%2d%2d')

%% load AARTFAAC covariane matrix
acc = acm + acm' - diag(diag(acm));
Nelem = size(acc, 1);

% load ITRF positions of antennas
fid = fopen('LBA_OUTER_AARTFAAC_config.txt', 'r');
posfiledata = textscan(fid, '%s %s [%f,%f,%f]');
posITRF = zeros(Nelem, 3);
posITRF(:, 1) = posfiledata{3};
posITRF(:, 2) = posfiledata{4};
posITRF(:, 3) = posfiledata{5};

% convert positions to local horizon frame
% rotation matrix taken from AntennaField.conf file from CS002
rotmat = [-0.1195950000, -0.7919540000, 0.5987530000; ...
           0.9928230000, -0.0954190000, 0.0720990000; ...
           0.0000330000,  0.6030780000, 0.7976820000];
poslocal = posITRF * rotmat;

% provide meta-data
l = -1:0.01:1;                            % (l,m)-grid for sky images
m = -1:0.01:1;
freq = 59672546;                           % frequency in Hz
%tobs = datenum([2011, 9, 21, f_tobs(1), f_tobs(2), f_tobs(3)]);  % time in UTC
tobs = datenum([2011, 9, 21, 13, 39, 8]); 

normal = [0.598753, 0.072099, 0.797682].'; % normal to CS002 field (ITRF)

% specify calibration sources from 3CR catalog (Cas A, Cyg A, Tau A, Vir A and Sun)
srcsel =  [324, 283, 88, 179, 0];

%% calibrate the data
[cal, sigmas, Sigman] = statcal(acc, tobs, freq, posITRF, srcsel, normal, 4, 20, eye(Nelem));
acccal = (cal' * cal) .* acc;

%% make sky maps
skymap = acm2skyimage(acc, poslocal(:, 1), poslocal(:, 2), freq, l, m);
skymapcal = acm2skyimage(acccal, poslocal(:, 1), poslocal(:, 2), freq, l, m);

%% show the results
mask = NaN(length(l));
mask(meshgrid(l).^2 + meshgrid(m).'.^2 < 1) = 1;
figure
imagesc(l, m, skymap .* mask);
set(gca, 'FontSize', 16);
axis equal
axis tight
xlabel('South \leftarrow m \rightarrow North');
ylabel('East \leftarrow l \rightarrow West');

figure
imagesc(l, m, skymapcal .* mask);
set(gca, 'FontSize', 16);
axis equal
axis tight
xlabel('South \leftarrow m \rightarrow North');
ylabel('East \leftarrow l \rightarrow West');

%% Save relevant variables
outfilename = strrep (tok, '.mat', '_var.mat');
save (outfilename, 'cal', 'sigmas', 'Sigman', 'skymap', 'skymapcal', 'tobs', 'freq', 'l', 'm');
