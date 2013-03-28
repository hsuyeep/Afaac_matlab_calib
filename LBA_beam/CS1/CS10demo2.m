close all
fig1h = figure;
fig2h = figure;

while (1)
dos('rsync -av lofarsys@129.125.99.51:sjw/xc-demo-20090403/*.dat ./demodata');

clear datafiles Mx My dirtymapx dirtymapy data

% parameters
c = 2.99792e8;
freq = 0:80e6/512:80e6-1e-3;
load CS1/CS10_config20080220.mat
dist = sqrt((meshgrid(xpos) - meshgrid(xpos).').^2 + (meshgrid(ypos) - meshgrid(ypos).').^2);
D = max(dist(:));
lat = 52 + 54/60 + 42.6/3600; % geographical location of CS010
lon = 6 + 52/60 + 3.17/3600;

% antenna selection for x-dipoles (from calibration measurement)
%sel = [1:18, 20:32 34:36 38 39 42:48];
selx = 1:48; %[1:14 16:36 38 39 42:48];

sely = 1:48; %[1:39 42:48];

% find data files
datadir = './demodata/';
[status, lsout] = dos(['ls -1 ' datadir '*.dat']);
datafiles = textscan(lsout, '%s\n');
nfiles = length(datafiles{1});
offset = length(datadir);

% check whether data is available
if (strcmp(datafiles{1}{1}(1:2), 'ls'))
    disp('No data available');
    pause(10);
    continue
end

% determine subbands
if (nfiles > 25)
    filesel = nfiles-24:nfiles;
else
    filesel = 1:nfiles;
end
sb = zeros(min(25, nfiles), 1);
for idx = 1:min(25, nfiles)
    sb(idx) = str2double(datafiles{1}{filesel(idx)}((19:21) + offset)) + 1;
end

% read and process data
% define grid based on highest available frequency
lambda = c / freq(max(sb));
dl = 0.75 * lambda / D;
l = 0:dl:1;
l = [-fliplr(l(2:end)), l];
[lgrid, mgrid] = meshgrid(l);
dist = sqrt(lgrid.^2 + mgrid.^2);
lgrid(dist > 1) = NaN;
mgrid(dist > 1) = NaN;
lgrid = lgrid(~isnan(lgrid));
mgrid = mgrid(~isnan(mgrid));
% initialize deconvolution matrix and dirty map vector
Mx = 0;
My = 0;
dirtymapx = 0;
dirtymapy = 0;
skymapx = 0;
skymapy = 0;
for idx = 1:length(sb);
    disp(sb(idx));
    
    % read data
    fid = fopen(datafiles{1}{filesel(idx)}, 'r');
    data = fread(fid, 'double');
    fclose(fid);
    R = reshape(data(1:2:end) + i * data(2:2:end), [96 96]);
    Rxx = R(1:2:end, 1:2:end);
    Rxx = Rxx(selx, selx);
    Ryy = R(2:2:end, 2:2:end);
    Ryy = Ryy(sely, sely);
    xposx = xpos(selx);
    yposx = ypos(selx);
    xposy = xpos(sely);
    yposy = ypos(sely);
    t_obs = datenum(datafiles{1}{filesel(idx)}(offset+1:offset+15), 'yyyymmdd_HHMMss');
    
    % FT imaging
    skymapx = skymapx + acm2skyimage(Rxx, xposx, yposx, freq(sb(idx)), l, l);
    skymapy = skymapy + acm2skyimage(Ryy, xposy, yposy, freq(sb(idx)), l, l);
    % calibrate data
    [calx, sigmasx, Sigmanx, lsrcx, msrcx] = statcal(Rxx, t_obs, freq(sb(idx)), xposx, yposx, lon, lat, 4);
    Rxx = (calx' * calx) .* Rxx;
    [caly, sigmasy, Sigmany, lsrcy, msrcy] = statcal(Ryy, t_obs, freq(sb(idx)), xposy, yposy, lon, lat, 4);
    Ryy = (caly' * caly) .* Ryy;
    
    % for Amir
    acc(:, :, idx) = Rxx;
    freqacc(idx) = freq(sb(idx));
    
    % compute dirty maps and deconvolution matrices for x- and y-arrays
    lambdasb = c / freq(sb(idx));
    A = exp(-(2 * pi * i / lambdasb) * (xposx * lgrid.' + yposx * mgrid.'));
    invR = inv(Rxx);
    Mx = Mx + abs(A' * invR * A).^2 - abs(A' * invR).^2 * inv(abs(invR).^2) * abs(invR * A).^2;
    dirtymapx = dirtymapx + real((conj(A' * invR) .* A') * ones(length(xposx), 1) - abs(A' * invR).^2 * inv(abs(invR).^2) * diag(invR));
    A = exp(-(2 * pi * i / lambdasb) * (xposy * lgrid.' + yposy * mgrid.'));
    invR = inv(Ryy);
    My = My + abs(A' * invR * A).^2 - abs(A' * invR).^2 * inv(abs(invR).^2) * abs(invR * A).^2;
    dirtymapy = dirtymapy + real((conj(A' * invR) .* A') * ones(length(xposy), 1) - abs(A' * invR).^2 * inv(abs(invR).^2) * diag(invR));
end

% deconvolve dirty map and show the result
figure(fig1h);
imagesc(l, l, flipud(rot90(skymapx + skymapy)))
set(gca, 'FontSize', 16, 'XDir', 'reverse', 'YDir', 'normal');
axis equal
xlabel('South');
ylabel('East');
set(gcf, 'Position', [3840 400 840 800]);
%set(gcf, 'Position', [0 400 840 800]);
figure(fig2h);
imdeconvx = zeros(size(dist));
imdeconvy = zeros(size(dist));
imdeconvx(dist < 1) = Mx \ dirtymapx;
imdeconvy(dist < 1) = My \ dirtymapy;
Imap = conv2(imdeconvx + imdeconvy, [0.25 0.5 0.25; 0.5 1 0.5; 0.25 0.5 0.25]);
imagesc(l, l, Imap);
caxis([min(Imap(:)), 0.8 * max(Imap(:))]);
set(gca, 'FontSize', 16, 'XDir', 'reverse', 'YDir', 'normal');
axis equal
xlabel('South');
ylabel('East');
set(gcf, 'Position', [4680 400 840 800]);
%set(gcf, 'Position', [840 400 840 800]);
drawnow
pause(0.5);

end
