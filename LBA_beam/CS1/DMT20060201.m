clear all
close all

alpha = (0:1:360) * (pi / 180);
delta = [-10 10 30 50 70 90] * (pi / 180);
t_obs = datenum(2006, 1, 1, 0, 0, 0);
lambda = 2.99792e8 / 240e6;
l0 = 0:0.05:1;
m0 = zeros(size(l0));
dl = 0.02;
lgrid = -1:dl:1;
mgrid = -1:dl:1;
spacing = 1.2;
nelem = 7;
gain = ones(nelem^2, 1);
pos = (-(nelem-1)/2:(nelem-1)/2) * spacing;
[xpos0, ypos0] = meshgrid(pos);
xpos = xpos0 * cos(pi/4) + ypos0 * sin(pi/4);
ypos = xpos0 * sin(pi/4) - ypos0 * cos(pi/4);

dist = sqrt((meshgrid(xpos(:)) - meshgrid(xpos(:)).').^2 + (meshgrid(ypos(:)) - meshgrid(ypos(:)).').^2);
M = 0.3 * 1.2 * dist.^-1 .* exp(2 * pi * i * dist / lambda);
M(isnan(M)) = gain + 0.0 * (randn(nelem^2, 1) + i * randn(nelem^2, 1)) / sqrt(2);
gain = M; %eye(nelem^2);
%load positions
%zooi;
%xpos = rand(nelem) * nelem * spacing - nelem * spacing / 2;
%ypos = rand(nelem) * nelem * spacing - nelem * spacing / 2;

%[tilex, tiley] = meshgrid(pos(1:4:end) + 1.5 * spacing);
%taperx = 1 - (abs(tilex) / (max(tilex(:)) + spacing)).^2;
%tapery = 1 - (abs(tiley) / (max(tiley(:)) + spacing)).^2;
%radius = sqrt(xpos(:).^2 + ypos(:).^2);
%gain = 1 - (radius / 45).^2;
%gain = kron(taperx .* tapery, ones(4));
%gain = taperx .* tapery;
%gain(1) = 0;

for idx = 1:length(delta)
    [l, m] = radectolm(alpha, delta(idx), JulianDay(t_obs));
    trackl(idx, :) = l;
    trackm(idx, :) = m;
end

df = 0e6;
fc = 2.99792e8 / lambda;
freq = fc-df/2:0.5e6:fc+df/2;
nfreq = length(freq);
phi = 0; %:pi/16:pi-1e-3;
ccbeam = zeros(length(l0), length(phi)^2, length(lgrid), length(mgrid));

fidx = 1;
for l0idx = 1:length(l0); %fidx = 1:nfreq;
disp(num2str(fidx));
pause(0.01);

for idx = 1:length(phi);
    tic
    disp(num2str(phi(idx) * 180 / pi));
    pause(0.01);
    arrayvbeam(idx, :, :) = xytolm(xpos(:) * sin(phi(idx)) + ypos(:) * cos(phi(idx)), xpos(:) * cos(phi(idx)) - ypos(:) * sin(phi(idx)), gain, lgrid, mgrid, 2.99792e8 / freq(fidx), l0(l0idx), m0(l0idx));
    toc
end
for idx1 = 1:length(phi)
    for idx2 = 1:length(phi)
        if (idx2 >= idx1)
            ccbeam(l0idx, length(phi) * (idx1-1) + idx2, :, :) = arrayvbeam(idx1, :, :) .* conj(arrayvbeam(idx2, :, :));
        end
        if (idx2 < idx1)
            ccbeam(l0idx, length(phi) * (idx1-1) + idx2, :, :) = conj(arrayvbeam(idx1, :, :)) .* arrayvbeam(idx2, :, :);
        end
    end
end
end % fidx=1:nfreq

save ccorbeam ccbeam
arraypbeam = real(squeeze(mean(mean(ccbeam, 1), 2)) / nelem^2);

dipolepbeam = calculateBeam(lgrid, mgrid, 30e6);
%elempbeam = ones(length(lgrid), length(mgrid));
elempbeam = dipolepbeam + dipolepbeam.';
dist = sqrt(meshgrid(lgrid).^2 + meshgrid(mgrid).'.^2);
mask = 1.0 * (dist < 1);
maskNaN = mask;
maskNaN(mask == 0) = NaN;

att = ones(length(0:dl:1), 1);
%att = delaysmearing(0:dl:1, df, 500);
smearing = interp1(0:dl:1-dl, att(1:end-1), dist);

figure
peakval = max(max(elempbeam .* mask));
elempbeam = elempbeam / peakval;
peakval = max(abs(arrayvbeam(:)));
arrayvbeam = arrayvbeam / peakval;
imagesc(lgrid, mgrid, 10 * log10(arraypbeam .* elempbeam .* smearing .* maskNaN));
set(gca, 'FontSize', 16);
caxis([-100, 0]);
xlabel('l');
ylabel('m');
set(colorbar, 'FontSize', 16);
%title('celestial tracks over statian beam pattern in dB')
hold on
h = plot(trackl.', trackm.');
set(h, 'Color', [1 1 1], 'LineWidth', 2);
axis equal
hold off

% main lobe sensitivity in Jy
sens0 = 50e-3 * (30e6/fc)^0.7;
limflux = (sens0 ./ (arraypbeam .* elempbeam .* smearing)) .* maskNaN;
srccount = zeros(size(limflux));
srccount(limflux<1e-3) = 53*limflux(limflux<1e-3).^(-1.2);
srccount((limflux>=1e-3)&(limflux<1e-1)) = 2023*limflux((limflux>=1e-3)&(limflux<1e-1)).^(-0.7);
srccount(limflux>1e-1) = 150*limflux(limflux>=1e-1).^(-1.5);

srperpixel = (dl^2 ./ sqrt(1 - dist.^2)) .* mask;
srcperpixel = srccount .* srperpixel;
srcperpixel(isnan(srcperpixel)) = 0;
senspat = arraypbeam .* elempbeam .* smearing .* maskNaN;
detected = zeros(size(limflux));
detected(senspat>0.1) = sum(srcperpixel(senspat>0.1));
detected((senspat>0.01)&(senspat<=0.1)) = sum(srcperpixel((senspat>0.01)&(senspat<=0.1)));
detected((senspat>0.001)&(senspat<=0.01)) = sum(srcperpixel((senspat>0.001)&(senspat<=0.01)));
detected((senspat>0.0001)&(senspat<=0.001)) = sum(srcperpixel((senspat>0.0001)&(senspat<=0.001)));
detected(senspat<=0.0001) = sum(srcperpixel(senspat<=0.0001));
%imagesc(lgrid, mgrid, log10(detected));
%set(gca, 'FontSize', 10);
%xlabel('l');
%ylabel('m');
%set(colorbar, 'FontSize', 10);

save result.mat