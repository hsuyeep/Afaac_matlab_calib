% Demonstration of AARTFAAC imaging
%
% SJW, 11 January 2012

%% start with a clean workspace
clear all
close all

% load AARTFAAC covariane matrix
load acm_sterp
% the array covariance matrix is stored in upper triangular form, so we
% complete it
acc = acm + acm' - diag(diag(acm));
Nelem = size(acc, 1);

% load ITRF positions of antennas
fid = fopen('LBA_OUTER_AARTFAAC_config.txt', 'r');
posfiledata = textscan(fid, '%s %s [%f,%f,%f]');
posITRF = zeros(Nelem, 3);
posITRF(:, 1) = posfiledata{3};
posITRF(:, 2) = posfiledata{4};
posITRF(:, 3) = posfiledata{5};
fclose(fid);

% convert positions to local horizon frame
% rotation matrix taken from AntennaField.conf file from CS002
rotmat = [-0.1195950000, -0.7919540000, 0.5987530000; ...
           0.9928230000, -0.0954190000, 0.0720990000; ...
           0.0000330000,  0.6030780000, 0.7976820000];
poslocal = posITRF * rotmat;

% provide meta-data
l = -1:0.005:1;                            % (l,m)-grid for sky images
m = -1:0.005:1;
c = 2.99792e8;                             % spped of light in m/s
lon = 6.869837540;                         % longitude of CS002 in degrees
lat = 52.915122495;                        % latitude of CS002 in degrees
freq = 59672546;                           % frequency in Hz
tobs = datenum([2011, 9, 21, 12, 39, 8]);  % time in UTC
normal = [0.598753, 0.072099, 0.797682].'; % normal to CS002 field (ITRF)

% specify calibration sources (Cas A, Cyg A, Tau A, Vir A and Sun)
srcsel =  [324, 283, 88, 179, 0];

%% calibrate the data
% whitening of the array covariance matrix (required for DOA estimation)
tic
acc = acc ./ sqrt(diag(acc) * diag(acc).');
% calibrate using source positions from 3C catalog and a baseline
% restriction of 10 wavelengths or 60 meters (whichever is smaller)
[cal, sigmas, Sigman] = statcal(acc, tobs, freq, posITRF, srcsel, normal, 10, 60, eye(Nelem));
disp('First round of antenna calibration completed');
toc

% Actual source positions may differ from the known positions due to
% ionospheric image distortions, so we use a high-resolution DOA estimation
% technique (weighted subspace fitting, WSF) to estimate source positions
% relative to the position of Cas A.
% We only estimate source positions for sources with an apparent flux
% larger than 1% of the apparent flux of Cas A.
tic
sel = sigmas > 0.01;
nsrc = sum(sel);
% use source positions from 3C catalog as initial estimate
load srclist3CR
if (sum(srcsel == 0) ~= 0)
    [raSun, decSun] = SunRaDec(JulianDay(tobs));                % J2000
    rasrc = [srclist3CR(srcsel(srcsel ~= 0)).alpha, raSun].';   % B1950
    decsrc = [srclist3CR(srcsel(srcsel ~= 0)).delta, decSun].'; % B1950
    epoch = [true(length(srcsel(srcsel ~= 0)), 1); false];
else
    rasrc = [srclist3CR(srcsel(srcsel ~= 0)).alpha];  % B1950
    decsrc = [srclist3CR(srcsel(srcsel ~= 0)).delta]; % B1950
    epoch = true(length(srcsel), 1);
end
srcpos0 = radectoITRF(rasrc(sel), decsrc(sel), epoch, JulianDay(tobs));
thsrc0 = asin(srcpos0(:, 3));
phisrc0 = atan2(srcpos0(:, 2), srcpos0(:, 1));

% find source positions using WSF algorithm
[v, d] = eig(acc);
[d, order] = sort(diag(abs(d)), 'descend');
v = v(:, order);
Wopt = (diag(d(1:nsrc)) - mean(diag(squeeze(Sigman))) * eye(nsrc))^2 / diag(d(1:nsrc));
Es = v(:, 1:nsrc);
EsWEs = Es * Wopt * Es';
theta = fminsearch(@(x) WSFcostfunITRF(x, EsWEs, diag(conj(1 ./ cal)), freq, posITRF), [phisrc0; thsrc0]);
phisrchat = theta(1:nsrc);
thsrchat = theta(nsrc+1:end);
srcposhat = [cos(phisrchat) .* cos(thsrchat), sin(phisrchat) .* cos(thsrchat), sin(thsrchat)];
disp('Position estimation using WSF done');
toc

% use updated source positions to improve calibration
tic
A = exp(-(2 * pi * 1i * freq / c) * (posITRF * srcposhat.'));
[ghat, sigmas2, Sigman2] = cal_ext(acc, A, sigmas(sel), squeeze(abs(Sigman)) > 0);
cal2 = (1 ./ ghat)';
% calibrate visibilities removing extended emission and system noise
acccal = (cal2' * cal2) .* (acc - squeeze(Sigman2));
disp('Calibration factors updated');
toc

%% A-team subtraction
% We start by removing the A-team sources. Since the Sun is a barely
% resolved source given the AARTFAAC resolution of 0.8 degrees at 60 MHz,
% it requires separate treatment.
% The difference between the actual positions as estimated by WSF and the
% positions in the catalog is sufficiently large to make subtraction based
% on catalog positions fail, we use the source positions estimated above.
% reduce srcsel to source that appear brighter than 1% if Cas A, note that
% 0 is used to denote the Sun
tic
srcsel = srcsel(sel);
% construct model visibilities for A-team sources
A = exp(-(2 * pi * 1i * freq / c) * (posITRF * srcposhat(srcsel ~= 0, :).'));
RAteam = A * diag(sigmas2(srcsel ~= 0)) * A';
% subtract A-team
accsubAteam = acccal - RAteam;
disp('A-team subtracted');
toc

%% Sun subtraction
% The Sun is modeled using a 3-sparse grid of pixels around the estimated
% position of the Sun. High resolution imaging using the MVDR beamformer to
% obtain a model of the Sun in the image domain proved unsuccessful.
% Finding the 3-sparse solution using an exhaustive search it the most
% computationally demanding part of this script. Hopefully Arash can
% improve on this using convex optimization.
tic
phi0sun = atan2(srcposhat(srcsel == 0, 2), srcposhat(srcsel == 0, 1));
th0sun = asin(srcposhat(srcsel == 0, 3));
% define grid of pixels
delta = 0.01;
res = 2 * delta / 8;
[phisun, thsun] = meshgrid(phi0sun-delta:res:phi0sun+delta, th0sun-delta:res:th0sun+delta);
phith_dist = sqrt((phisun(:) - phi0sun).^2 + (thsun(:) - th0sun).^2);
phisun = phisun(phith_dist < delta);
thsun = thsun(phith_dist < delta);
sunpos = [cos(phisun) .* cos(thsun), sin(phisun) .* cos(thsun), sin(thsun)];

% solar component estimation using sparse reconstruction
A = exp(-(2 * pi * 1i * freq / c) * (posITRF * sunpos.'));
fluxSun = 1.20 * sigmas2(srcsel == 0);
cvx_begin
    variable sigmaSun(size(A, 2));
    minimize norm(accsubAteam(:) - khatrirao(conj(A), A) * sigmaSun, 2);
    subject to
    norm(sigmaSun, 1) <= fluxSun;
cvx_end

% iterate, increasing resolution
sigmaSun(sigmaSun < 1e-3) = 0;
phi0sun = phisun(sigmaSun ~= 0);
th0sun = thsun(sigmaSun ~= 0);
res2 = res / 3;
[phigrid, thgrid] = meshgrid(-res2:res2:res2);
phisun2 = kron(phi0sun, ones(9, 1)) + kron(ones(length(phi0sun), 1), phigrid(:));
thsun2 = kron(th0sun, ones(9, 1)) + kron(ones(length(th0sun), 1), thgrid(:));
sunpos2 = [cos(phisun2) .* cos(thsun2), sin(phisun2) .* cos(thsun2), sin(thsun2)];
A = exp(-(2 * pi * 1i * freq / c) * (posITRF * sunpos2.'));
fluxSun = 1.20 * sigmas2(srcsel == 0);
cvx_begin
    variable sigmaSun(size(A, 2));
    minimize norm(accsubAteam(:) - khatrirao(conj(A), A) * sigmaSun, 2);
    subject to
    norm(sigmaSun, 1) <= fluxSun;
cvx_end

% try subtraction
sigmaSun(sigmaSun < 1e-3) = 0;
RSunS = A * diag(sigmaSun) * A';
accsubSunS = acccal - RAteam - RSunS;
disp('Sun subtracted');
toc

%% initiallize matrix containing least squares cost value
%cost = NaN(length(lsun(:)), length(lsun(:)), length(lsun(:)));
% pre-calculation, system noise added to ensure the covariance matrix is PD
%invR = inv(accsubAteam + (cal2' * cal2) .* diag(diag(Sigman2)));
% Perform exhaustive search over all possible combinations of three
% components on the grid defined above. Deconvolution is done using LS
% imaging.
%for idx1 = 1:length(lsun(:))-2
%    disp(['current value of idx1: ' num2str(idx1)]);
%    pause(0.01); % to ensure progress indication is displayed
%    for idx2 = idx1+1:length(lsun(:))-1
%        for idx3 = idx2+1:length(lsun(:))
%            A = exp(-(2 * pi * 1i * freq / c) * (poslocal(:, 1) * [lsun(idx1), lsun(idx2), lsun(idx3)] + poslocal(:, 2) * [msun(idx1), msun(idx2), msun(idx3)]));
%            % calculate deconvolution matrix
%            M = abs(A' * invR * A).^2;
%            % calculate dirty image
%            sd = khatrirao(conj(invR * A), invR * A)' * accsubAteam(:);
%            % estimate flux
%            flux = real(M \ sd);
%            % if flux values are physically meaningful, assess quality of
%            % fit using least squares cost function
%            if sum(flux < 0) == 0
%                cost(idx1, idx2, idx3) = sum(sum(abs(accsubAteam - A * diag(flux) * A').^2));
%            end
%        end
%    end
%end
% find the three components giving lowest least squares cost
%[~, minidx] = min(cost(:));
%[idx1, idx2, idx3] = ind2sub([length(lsun(:)), length(lsun(:)), length(lsun(:))], minidx);
% construct source model using these compoents
%A = exp(-(2 * pi * 1i * freq / c) * (poslocal(:, 1) * [lsun(idx1), lsun(idx2), lsun(idx3)] + poslocal(:, 2) * [msun(idx1), msun(idx2), msun(idx3)]));
%M = abs(A' * invR * A).^2;
%sd = khatrirao(conj(invR * A), invR * A)' * accsubAteam(:);
%flux = M \ sd;
%RSun = A * diag(flux) * A';
% subtract A-team and Sun
%accsubSun = acccal - RAteam - RSun;

%% all sky image after calibration, but before source subtraction
u = meshgrid(poslocal(:, 1)) - meshgrid(poslocal(:, 1)).';
v = meshgrid(poslocal(:, 2)) - meshgrid(poslocal(:, 2)).';
duv = 1;
Nuv = 1000;
uvsize = 1024;
[skymapcal, ~] = fft_imager(acccal + (cal2' * cal2) .* squeeze(Sigman2), u, v, duv, Nuv, uvsize);
lambda = c / freq;
ltick_res = lambda / (uvsize * duv);
ltick_max = uvsize * ltick_res / 2;
ltick = -ltick_max:ltick_res:ltick_max-1e-3;
%skymapcal = acm2skyimage(acccal + (cal2' * cal2) .* squeeze(Sigman2), poslocal(:, 1), poslocal(:, 2), freq, l, m);
mask = NaN(length(ltick));
mask(meshgrid(ltick).^2 + meshgrid(ltick).'.^2 < 1) = 1;
figure
imagesc(ltick, ltick, real(skymapcal) .* mask);
set(gca, 'FontSize', 16);
title('calibrated sky map, extended emission removed');
axis equal
axis tight
xlabel('South \leftarrow m \rightarrow North');
ylabel('East \leftarrow xl \rightarrow West');
set(colorbar, 'FontSize', 16);
axis([-1 1 -1 1]);

%% all-sky map after removal of the A-team sources
[skymapsubAteam, ~] = fft_imager(accsubAteam, u, v, duv, Nuv, uvsize);
%skymapsubAteam = acm2skyimage(accsubAteam, poslocal(:, 1), poslocal(:, 2), freq, l, m);
figure
imagesc(ltick, ltick, real(skymapsubAteam) .* mask);
set(gca, 'FontSize', 16);
title('extended emission and A-team removed');
axis equal
axis tight
xlabel('South \leftarrow m \rightarrow North');
ylabel('East \leftarrow l \rightarrow West');
set(colorbar, 'FontSize', 16);
axis([-1 1 -1 1]);

%% all-sky map after removal of the A-team sources and the Sun
%mask = NaN(length(l));
%mask(meshgrid(l).^2 + meshgrid(m).'.^2 < 1) = 1;
%uvdist = sqrt(u.^2 + v.^2);
%uvweight = min(uvdist / 60, 1);
[skymapsubSun, ~] = fft_imager(accsubSunS, u, v, duv, Nuv, uvsize);
%skymapsubSun = acm2skyimage(accsubSunS, poslocal(:, 1), poslocal(:, 2), freq, l, m);
figure
imagesc(ltick, ltick, real(skymapsubSun) .* mask);
set(gca, 'FontSize', 16);
title('extended emission, Sun and A-team removed');
axis equal
axis tight
xlabel('South \leftarrow m \rightarrow North');
ylabel('East \leftarrow l \rightarrow West');
set(colorbar, 'FontSize', 16);
axis([-1 1 -1 1]);
caxis([-500 5500])