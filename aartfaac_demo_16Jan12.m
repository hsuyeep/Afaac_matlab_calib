% Demonstration of AARTFAAC imaging
%
% SJW, 11 January 2012

%% start with a clean workspace
clear all
close all

% load AARTFAAC covariane matrix
load acm_sterp
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
l = -1:0.005:1;                            % (l,m)-grid for sky images
m = -1:0.005:1;
freq = 59672546;                           % frequency in Hz
tobs = datenum([2011, 9, 21, 12, 39, 8]);  % time in UTC
normal = [0.598753, 0.072099, 0.797682].'; % normal to CS002 field (ITRF)

% specify calibration sources (Cas A, Cyg A, Tau A, Vir A and Sun)
srcsel =  [324, 283, 88, 179, 0];

%% calibrate the data
[cal, sigmas, Sigman] = statcal(acc, tobs, freq, posITRF, srcsel, normal, 10, 60, eye(Nelem));
acccal = (cal' * cal) .* (acc - squeeze(Sigman));

%% model of bright sources
c = 2.99792e8;
load srclist3CR
if (sum(srcsel == 0) ~= 0)
    [raSun, decSun] = SunRaDec(JulianDay(tobs));                % J2000
    rasrc = [srclist3CR(srcsel(srcsel ~= 0)).alpha, raSun].';   % B1950
    decsrc = [srclist3CR(srcsel(srcsel ~= 0)).delta, decSun].'; % B1950
    epoch = [true(length(srcsel(srcsel ~= 0)), 1); false];
else
    rasrc = srclist3CR(srcsel).alpha;  % B1950
    decsrc = srclist3CR(srcsel).delta; % B1950
    epoch = true(length(srcsel), 1);
end
srcpos = radectoITRF(rasrc, decsrc, epoch, JulianDay(tobs));
up = srcpos * normal > 0;
A = exp(-(2 * pi * 1i * freq / c) * (posITRF * srcpos(up, :).'));
Sigma = diag(sigmas);
Rsrc = A * Sigma * A';
accsubtract = acccal - Rsrc;

%% subtract A-team sources using estimated positions
c = 2.99792e8;
lsrchat = [0.2225; 0.736; -0.15];
msrchat = [0.888; 0.534; -0.6435];
up = [true true false true];
sigmas0 = sigmas(1:4);
Sigma = diag(sigmas0(up));
A = exp(-(2 * pi * 1i * freq / c) * (poslocal(:, 1) * lsrchat.' + poslocal(:, 2) * msrchat.'));
RAteam = A * Sigma * A';
accsubtract2 = acccal - RAteam;

%% separate model for the Sun using projected baselines
l0sun = -0.3102;
m0sun = -0.752;
A = exp(-(2 * pi * 1i * freq / c) * (poslocal(:, 1) * l0sun + poslocal(:, 2) * m0sun));
u = meshgrid(poslocal(:, 1)) - meshgrid(poslocal(:, 1)).';
v = meshgrid(poslocal(:, 2)) - meshgrid(poslocal(:, 2)).';
uproj = u * sqrt(1 - l0sun^2);
vproj = v * sqrt(1 - m0sun^2);
projbaseline = sqrt(uproj.^2 + vproj.^2);
RSun = 1.39 * (A * sigmas(5) * A') .* exp(-projbaseline / (1.5 * max(projbaseline(:))));
accsubtract3 = acccal - RSun - RAteam;
skymapnew3 = acm2skyimage(accsubtract3, poslocal(:, 1), poslocal(:, 2), freq, l, m);
imagesc(l, m, skymapnew3); set(colorbar, 'FontSize', 16);

%% op hoge resolutie wegprojecteren van de zon
l0sun = -0.3115;
m0sun = -0.751;
lpgrid = l0sun-0.01:1e-4:l0sun+0.01;
mpgrid = m0sun-0.01:1e-4:m0sun+0.01;
Aproj = exp(-(2 * pi * 1i * freq / c) * (poslocal(:, 1) * lpgrid + poslocal(:, 2) * mpgrid));
[eigvec, eigval] = eig(Aproj' * Aproj);
invidx = (diag(eigval) > max(diag(eigval)) / 1e4);
invmat = eigvec(:, invidx) * diag(1 ./ diag(eigval(invidx, invidx))) * eigvec(:, invidx)';
Pproj = eye(size(Aproj, 1)) - Aproj * (invmat * Aproj');
accproj = Pproj * (acccal - RAteam) * Pproj';
skymapp = acm2skyimage(accproj, poslocal(:, 1), poslocal(:, 2), freq, l, m);

%% make sky maps
%skymap = acm2skyimage(acc, poslocal(:, 1), poslocal(:, 2), freq, l, m);
skymapcal = acm2skyimage(acccal, poslocal(:, 1), poslocal(:, 2), freq, l, m);

%% show the results
mask = NaN(length(l));
mask(meshgrid(l).^2 + meshgrid(m).'.^2 < 1) = 1;
figure
imagesc(l, m, skymapcal .* mask);
set(gca, 'FontSize', 16);
title('calibrated sky map');
axis equal
axis tight
xlabel('South \leftarrow m \rightarrow North');
ylabel('East \leftarrow l \rightarrow West');

figure
imagesc(l, m, skymapp .* mask);
set(gca, 'FontSize', 16);
title('A-team and extended emission subtracted, Sun projected out');
axis equal
axis tight
xlabel('South \leftarrow m \rightarrow North');
ylabel('East \leftarrow l \rightarrow West');
caxis([-500, 6000]);

figure
imagesc(l, m, skymapnew3 .* mask);
set(gca, 'FontSize', 16);
title('A-team, Sun and extended emission subtracted');
axis equal
axis tight
xlabel('South \leftarrow m \rightarrow North');
ylabel('East \leftarrow l \rightarrow West');
