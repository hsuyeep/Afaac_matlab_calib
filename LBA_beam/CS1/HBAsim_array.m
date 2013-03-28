clear all
close all

% parameter section
l = -1:0.01:1;
m = -1:0.01:1;
zenithangle = 40 * pi / 180;
azimuth = (0:6:179) * pi / 180;
l0 = 0.7;
m0 = 0.4;
lambda = 2.7:-0.025:1.2;
spacing = 1.25;
nelem1dim = 20;
pos = (-(nelem1dim-1)/2:(nelem1dim-1)/2) * spacing;
[xpos0, ypos0] = meshgrid(pos);
phi_array = pi * 45 / 180;
xpos = xpos0(:) * cos(phi_array) + ypos0(:) * sin(phi_array);
ypos = xpos0(:) * sin(phi_array) - ypos0(:) * cos(phi_array);
Mampl = 0.07;
Mdir = 0.9;
iter = 10;

% precomputations
% distances
distx = meshgrid(xpos) - meshgrid(xpos).';
disty = meshgrid(ypos) - meshgrid(ypos).';
dist = sqrt(distx.^2 + disty.^2);
% baseline orientations
phir = atan2(disty, distx);

% memory allocation
arraysignal = cell(length(lambda), length(azimuth));

for idx = 1:length(lambda)
    disp(['lambda = ' num2str(lambda(idx))]);
    
    % coupling matrix
    tic
    M = Mampl * (dist / lambda(idx)).^-1 .* (1 - Mdir * abs(cos(phir))) .* exp(2 * pi * i * dist / lambda(idx));
    M(~isfinite(M)) = 1;
    Meff = M;
    toc
    tic
    for pow = 1:iter
        Meff = Meff + (M - eye(length(xpos)))^pow;
    end
    toc
    tic
    for rotidx = 1:length(azimuth)
        % rotated elementpositions
        xpos = xpos0(:) * cos(phi_array + azimuth(rotidx)) + ypos0(:) * sin(phi_array + azimuth(rotidx));
        ypos = xpos0(:) * sin(phi_array + azimuth(rotidx)) - ypos0(:) * cos(phi_array + azimuth(rotidx));
        arraysignal{idx, rotidx} = squeeze(xytolm(xpos, ypos, Meff, l, m, lambda(idx), l0, m0));
    end
    toc
    tic
    save arraypattern
    toc
    disp([num2str(idx/length(lambda) * 100) '% complete']);
end
