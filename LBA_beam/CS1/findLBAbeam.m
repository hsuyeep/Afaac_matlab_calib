clear all
close all

lat = 52 + 54/60 + 42.6/3600; % geographical location of CS010
lon = 6 + 52/60 + 3.17/3600;

load srclist3C
srcsel = [324, 283]; % Cas A and Cyg A

load session20080309.mat
Nobs = size(sigmaxtot, 1);
Nsb = size(sigmaxtot, 2);

tobs = zeros(Nobs, Nsb);
lsrc = zeros(Nobs, Nsb, 2);
msrc = zeros(Nobs, Nsb, 2);
for idx = 1:Nobs
    disp(num2str(idx));
    eval(['load LBAmonitoring/snapshot' num2str(idx, '%03d'), '.mat']);
    eval(['tobs(idx, :) = data' num2str(idx, '%03d') '.t_obs;']);
    eval(['clear data' num2str(idx,'%03d')]);
    for sbidx = 1:512
        [lsrc(idx, sbidx, :), msrc(idx, sbidx, :)] = radectolm(srclist3C.alpha(srcsel), srclist3C.delta(srcsel), JulianDay(tobs(idx, sbidx)), lon, lat);
    end
end

thetasrc = asin(sqrt(lsrc.^2 + msrc.^2));
phisrc = atan2(lsrc, msrc);

sb = 250;
selobsx = sigmaxtot(:, sb, 2) > 0;
selobsy = sigmaytot(:, sb, 2) > 0;
phirot = pi/4;
orderth = 6;
orderphi = 3;
cthCasx = cos((0:(orderth-1)).' * squeeze(thetasrc(selobsx, sb, 1)).');
cthCygx = cos((0:(orderth-1)).' * squeeze(thetasrc(selobsx, sb, 2)).');
cphiCasx = 1 + cos(2 * (0:(orderphi-1)).' * (squeeze(phisrc(selobsx, sb, 1)).' + phirot));
sphiCasx = 1 + sin(2 * (1:(orderphi-1)).' * (squeeze(phisrc(selobsx, sb, 1)).' + phirot));
cphiCygx = 1 + cos(2 * (0:(orderphi-1)).' * (squeeze(phisrc(selobsx, sb, 2)).' + phirot));
sphiCygx = 1 + sin(2 * (1:(orderphi-1)).' * (squeeze(phisrc(selobsx, sb, 2)).' + phirot));

cthCasy = cos((0:(orderth-1)).' * squeeze(thetasrc(selobsy, sb, 1)).');
cthCygy = cos((0:(orderth-1)).' * squeeze(thetasrc(selobsy, sb, 2)).');
cphiCasy = 1 + cos(2 * (0:(orderphi-1)).' * (squeeze(phisrc(selobsy, sb, 1)).' - phirot));
sphiCasy = 1 + sin(2 * (1:(orderphi-1)).' * (squeeze(phisrc(selobsy, sb, 1)).' - phirot));
cphiCygy = 1 + cos(2 * (0:(orderphi-1)).' * (squeeze(phisrc(selobsy, sb, 2)).' - phirot));
sphiCygy = 1 + sin(2 * (1:(orderphi-1)).' * (squeeze(phisrc(selobsy, sb, 2)).' - phirot));

beamCas = [khatrirao(cthCasx, [cphiCasx; sphiCasx]), ...
           khatrirao(cthCasy, [cphiCasy; sphiCasy])];
beamCyg = [khatrirao(cthCygx, [cphiCygx; sphiCygx]), ...
           khatrirao(cthCygy, [cphiCygy; sphiCygy])];

[v, d] = eig(pinv(beamCyg.') * diag([sigmaxtot(selobsx, sb, 2); sigmaytot(selobsy, sb, 2)]) * beamCas.');

beamcoeff = v(:, 1) + v(:, 2);
[l, m] = meshgrid(-1:0.02:1, -1:0.02:1);
lmdist = sqrt(l.^2 + m.^2);
mask = ones(size(l));
mask(lmdist > 0.95) = NaN;
theta = asin(sqrt(l(:).^2 + m(:).^2));
phi = atan2(l(:), m(:));
ctheta = cos((0:(orderth-1)).' * theta.');
cphi = 1 + cos(2 * (0:(orderphi-1)).' * (phi.' + phirot));
sphi = 1 + sin(2 * (1:(orderphi-1)).' * (phi.' + phirot));
beam = beamcoeff' * khatrirao(ctheta, [cphi; sphi]);
imagesc(reshape(abs(beam), [101 101]) .* mask);
