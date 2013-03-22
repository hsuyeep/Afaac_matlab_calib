% Script to generate images from various stations, with ACMs already 
% assumed extracted via casapy and cleaned of external headers
% pep, Jan12
l = [-1:0.02:1]; m = [-1:0.02:1];

[posxpol_cs003, posypol_cs003, lon_cs003, lat_cs003] = parseAntennaField ('./AntennaFields/AntennaFieldCS003.conf', 1, 0);
CS003_avg_acc = fill_acm ('../cookdat/cs003_station_images/CS003_10s_60ch_acm_clean.txt', 47);
sky_CS003_avg = acm2skyimage (CS003_avg_acc, posxpol_cs003(:,1), posxpol_cs003(:,2), 59672546, l, m);

[posxpol_cs005, posypol_cs005, lon_cs005, lat_cs005] = parseAntennaField ('./AntennaFields/AntennaFieldCS005.conf', 1, 0);
CS005_avg_acc = fill_acm ('../cookdat/cs003_station_images/CS005_10s_60ch_acm_clean.txt', 47);
sky_CS005_avg = acm2skyimage (CS005_avg_acc, posxpol_cs005(:,1), posxpol_cs005(:,2), 59672546, l, m);

[posxpol_cs007, posypol_cs007, lon_cs007, lat_cs007] = parseAntennaField ('./AntennaFields/AntennaFieldCS007.conf', 1, 0);
CS007_avg_acc = fill_acm ('../cookdat/cs003_station_images/CS007_10s_60ch_acm_clean.txt', 47);
sky_CS007_avg = acm2skyimage (CS007_avg_acc, posxpol_cs007(:,1), posxpol_cs007(:,2), 59672546, l, m);

% show the result
dist = sqrt(meshgrid(l).^2 + meshgrid(m).'.^2);
mask = NaN(size(dist));
mask(dist < 1) = 1;
figure
set(gcf, 'Position', [0, 500, 1120 840]);
subplot(2, 2, 1);
imagesc(l, m, mean(sky_CS003_avg, 3).' .* mask);
set(gca, 'FontSize', 16, 'YDir', 'normal', 'XDir', 'reverse');
set(colorbar, 'FontSize', 16);
xlabel('East \leftarrow l \rightarrow West');
ylabel('South \leftarrow m \rightarrow North');
title('CS003, 10s integ, 60chans');

subplot(2, 2, 2);
imagesc(l, m, mean(sky_CS005_avg, 3).' .* mask);
set(gca, 'FontSize', 16, 'YDir', 'normal', 'XDir', 'reverse');
set(colorbar, 'FontSize', 16);
xlabel('East \leftarrow l \rightarrow West');
ylabel('South \leftarrow m \rightarrow North');
title('CS005, 10s integ, 60 chans');

subplot(2, 2, 3);
imagesc(l, m, mean(sky_CS007_avg, 3).' .* mask);
set(gca, 'FontSize', 16, 'YDir', 'normal', 'XDir', 'reverse');
set(colorbar, 'FontSize', 16);
xlabel('East \leftarrow l \rightarrow West');
ylabel('South \leftarrow m \rightarrow North');
title('CS007, 10s integ, 60 chans');
