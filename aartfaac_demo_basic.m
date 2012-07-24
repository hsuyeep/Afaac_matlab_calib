% Demonstration of AARTFAAC imaging
%
% SJW, 11 January 2012

%% start with a clean workspace
clear all
close all

% load AARTFAAC covariane matrix
load acm_sterp
acc = acm + acm' - diag(diag(acm)); % Remove autocorrelations
Nelem = size(acc, 1);

% load ITRF positions of antennas
%fid = fopen('LBA_OUTER_AARTFAAC_config.txt', 'r');
%posfiledata = textscan(fid, '%s %s [%f,%f,%f]');
%posITRF = zeros(Nelem, 3);
%posITRF(:, 1) = posfiledata{3};
%posITRF(:, 2) = posfiledata{4};
%posITRF(:, 3) = posfiledata{5};

% convert positions to local horizon frame
% rotation matrix taken from AntennaField.conf file from CS002
rotmat = [-0.1195950000, -0.7919540000, 0.5987530000; ...
           0.9928230000, -0.0954190000, 0.0720990000; ...
           0.0000330000,  0.6030780000, 0.7976820000];
%poslocal = posITRF * rotmat;
load poslocal.mat

% provide meta-data
l = -1:0.005:1;                            % (l,m)-grid for sky images
m = -1:0.005:1;
freq = 59672546;                           % frequency in Hz
tobs = datenum([2011, 9, 21, 12, 39, 8]);  % time in UTC
normal = [0.598753, 0.072099, 0.797682].'; % normal to CS002 field (ITRF)

% specify calibration sources (Cas A, Cyg A, Tau A, Vir A and Sun)
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

