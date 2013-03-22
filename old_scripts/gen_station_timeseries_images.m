% Script to generate images from various timeslices of a single stations, with ACMs already 
% assumed extracted vis casapy and cleaned of external headers
l = [-1:0.02:1]; m = [-1:0.02:1];

[posxpol_cs003, posypol_cs003, lon_cs003, lat_cs003] = parseAntennaField ('./AntennaFields/AntennaFieldCS003.conf', 1, 0);
CS003_ts0_avg_acc = fill_acm ('../cookdat/TSERIES/SB000_115900.txt', 47);
sky_CS003_ts0_avg = acm2skyimage (CS003_ts0_avg_acc, posxpol_cs003(:,1), posxpol_cs003(:,2), 59672546, l, m);

CS003_ts1_avg_acc = fill_acm ('../cookdat/TSERIES/SB000_120908.txt', 47);
sky_CS003_ts1_avg = acm2skyimage (CS003_ts1_avg_acc, posxpol_cs003(:,1), posxpol_cs003(:,2), 59672546, l, m);

CS003_ts2_avg_acc = fill_acm ('../cookdat/TSERIES/SB000_121908.txt', 47);
sky_CS003_ts2_avg = acm2skyimage (CS003_ts2_avg_acc, posxpol_cs003(:,1), posxpol_cs003(:,2), 59672546, l, m);

CS003_ts3_avg_acc = fill_acm ('../cookdat/TSERIES/SB000_122908.txt', 47);
sky_CS003_ts3_avg = acm2skyimage (CS003_ts3_avg_acc, posxpol_cs003(:,1), posxpol_cs003(:,2), 59672546, l, m);

CS003_ts4_avg_acc = fill_acm ('../cookdat/TSERIES/SB000_123908.txt', 47);
sky_CS003_ts4_avg = acm2skyimage (CS003_ts4_avg_acc, posxpol_cs003(:,1), posxpol_cs003(:,2), 59672546, l, m);

CS003_ts5_avg_acc = fill_acm ('../cookdat/TSERIES/SB000_124908.txt', 47);
sky_CS003_ts5_avg = acm2skyimage (CS003_ts5_avg_acc, posxpol_cs003(:,1), posxpol_cs003(:,2), 59672546, l, m);

CS003_ts6_avg_acc = fill_acm ('../cookdat/TSERIES/SB000_125908.txt', 47);
sky_CS003_ts6_avg = acm2skyimage (CS003_ts6_avg_acc, posxpol_cs003(:,1), posxpol_cs003(:,2), 59672546, l, m);

%CS003_ts7_avg_acc = fill_acm ('./SB000_130908.txt', 47);
%sky_CS003_ts7_avg = acm2skyimage (CS003_ts7_avg_acc, posxpol_cs003(:,1), posxpol_cs003(:,2), 59672546, l, m);

CS003_ts8_avg_acc = fill_acm ('../cookdat/TSERIES/SB000_131908.txt', 47);
sky_CS003_ts8_avg = acm2skyimage (CS003_ts8_avg_acc, posxpol_cs003(:,1), posxpol_cs003(:,2), 59672546, l, m);

%CS003_ts9_avg_acc = fill_acm ('./SB000_132908.txt', 47);
%sky_CS003_ts9_avg = acm2skyimage (CS003_ts9_avg_acc, posxpol_cs003(:,1), posxpol_cs003(:,2), 59672546, l, m);

CS003_ts10_avg_acc = fill_acm ('../cookdat/TSERIES/SB000_133908.txt', 47);
sky_CS003_ts10_avg = acm2skyimage (CS003_ts10_avg_acc, posxpol_cs003(:,1), posxpol_cs003(:,2), 59672546, l, m);

CS003_ts11_avg_acc = fill_acm ('../cookdat/TSERIES/SB000_134908.txt', 47);
sky_CS003_ts11_avg = acm2skyimage (CS003_ts11_avg_acc, posxpol_cs003(:,1), posxpol_cs003(:,2), 59672546, l, m);

CS003_ts12_avg_acc = fill_acm ('../cookdat/TSERIES/SB000_135908.txt', 47);
sky_CS003_ts12_avg = acm2skyimage (CS003_ts12_avg_acc, posxpol_cs003(:,1), posxpol_cs003(:,2), 59672546, l, m);

CS003_ts13_avg_acc = fill_acm ('../cookdat/TSERIES/SB000_140908.txt', 47);
sky_CS003_ts13_avg = acm2skyimage (CS003_ts13_avg_acc, posxpol_cs003(:,1), posxpol_cs003(:,2), 59672546, l, m);

CS003_ts14_avg_acc = fill_acm ('../cookdat/TSERIES/SB000_141908.txt', 47);
sky_CS003_ts14_avg = acm2skyimage (CS003_ts14_avg_acc, posxpol_cs003(:,1), posxpol_cs003(:,2), 59672546, l, m);

CS003_ts15_avg_acc = fill_acm ('../cookdat/TSERIES/SB000_142908.txt', 47);
sky_CS003_ts15_avg = acm2skyimage (CS003_ts15_avg_acc, posxpol_cs003(:,1), posxpol_cs003(:,2), 59672546, l, m);

CS003_ts16_avg_acc = fill_acm ('../cookdat/TSERIES/SB000_143908.txt', 47);
sky_CS003_ts16_avg = acm2skyimage (CS003_ts16_avg_acc, posxpol_cs003(:,1), posxpol_cs003(:,2), 59672546, l, m);

CS003_ts17_avg_acc = fill_acm ('../cookdat/TSERIES/SB000_144908.txt', 47);
sky_CS003_ts17_avg = acm2skyimage (CS003_ts17_avg_acc, posxpol_cs003(:,1), posxpol_cs003(:,2), 59672546, l, m);

CS003_ts18_avg_acc = fill_acm ('../cookdat/TSERIES/SB000_145908.txt', 47);
sky_CS003_ts18_avg = acm2skyimage (CS003_ts18_avg_acc, posxpol_cs003(:,1), posxpol_cs003(:,2), 59672546, l, m);

% show the result
dist = sqrt(meshgrid(l).^2 + meshgrid(m).'.^2);
mask = NaN(size(dist));
mask(dist < 1) = 1;
figure (1);
set(gcf, 'Position', [0, 500, 1120 840]);
%subplot(2, 2, 1);
imagesc(l, m, mean(sky_CS003_ts0_avg, 3).' .* mask);
set(gca, 'FontSize', 16, 'YDir', 'normal', 'XDir', 'reverse');
set(colorbar, 'FontSize', 16);
xlabel('East \leftarrow l \rightarrow West');
ylabel('South \leftarrow m \rightarrow North');
title('CS003, 10s integ, 40chans, 11:59');

%subplot(2, 2, 2);
figure(2);
imagesc(l, m, mean(sky_CS003_ts1_avg, 3).' .* mask);
set(gca, 'FontSize', 16, 'YDir', 'normal', 'XDir', 'reverse');
set(colorbar, 'FontSize', 16);
xlabel('East \leftarrow l \rightarrow West');
ylabel('South \leftarrow m \rightarrow North');
title('CS003, 10s integ, 40 chans, 12:09');

%subplot(2, 2, 3);
figure (3);
imagesc(l, m, mean(sky_CS003_ts2_avg, 3).' .* mask);
set(gca, 'FontSize', 16, 'YDir', 'normal', 'XDir', 'reverse');
set(colorbar, 'FontSize', 16);
xlabel('East \leftarrow l \rightarrow West');
ylabel('South \leftarrow m \rightarrow North');
title('CS003, 10s integ, 40chans, 12:19');

%subplot(2, 2, 4);
figure (4);
imagesc(l, m, mean(sky_CS003_ts3_avg, 3).' .* mask);
set(gca, 'FontSize', 16, 'YDir', 'normal', 'XDir', 'reverse');
set(colorbar, 'FontSize', 16);
xlabel('East \leftarrow l \rightarrow West');
ylabel('South \leftarrow m \rightarrow North');
title('CS003, 10s integ, 40 chans, 12:29');

%figure (2);
figure (5);
set(gcf, 'Position', [0, 500, 1120 840]);
%subplot(2, 2, 1);
imagesc(l, m, mean(sky_CS003_ts4_avg, 3).' .* mask);
set(gca, 'FontSize', 16, 'YDir', 'normal', 'XDir', 'reverse');
set(colorbar, 'FontSize', 16);
xlabel('East \leftarrow l \rightarrow West');
ylabel('South \leftarrow m \rightarrow North');
title('CS003, 10s integ, 40 chans, 12:39');

%subplot(2, 2, 2);
figure (6);
imagesc(l, m, mean(sky_CS003_ts5_avg, 3).' .* mask);
set(gca, 'FontSize', 16, 'YDir', 'normal', 'XDir', 'reverse');
set(colorbar, 'FontSize', 16);
xlabel('East \leftarrow l \rightarrow West');
ylabel('South \leftarrow m \rightarrow North');
title('CS003, 10s integ, 40 chans, 12:49');

%subplot(2, 2, 3);
figure (7);
imagesc(l, m, mean(sky_CS003_ts6_avg, 3).' .* mask);
set(gca, 'FontSize', 16, 'YDir', 'normal', 'XDir', 'reverse');
set(colorbar, 'FontSize', 16);
xlabel('East \leftarrow l \rightarrow West');
ylabel('South \leftarrow m \rightarrow North');
title('CS003, 10s integ, 40chans, 12:59');

%subplot(2, 2, 4);
figure (8);
imagesc(l, m, mean(sky_CS003_ts8_avg, 3).' .* mask);
set(gca, 'FontSize', 16, 'YDir', 'normal', 'XDir', 'reverse');
set(colorbar, 'FontSize', 16);
xlabel('East \leftarrow l \rightarrow West');
ylabel('South \leftarrow m \rightarrow North');
title('CS003, 10s integ, 40 chans, 13:19');



figure (9);
set(gcf, 'Position', [0, 500, 1120 840]);
%subplot(2, 2, 1);
imagesc(l, m, mean(sky_CS003_ts10_avg, 3).' .* mask);
set(gca, 'FontSize', 16, 'YDir', 'normal', 'XDir', 'reverse');
set(colorbar, 'FontSize', 16);
xlabel('East \leftarrow l \rightarrow West');
ylabel('South \leftarrow m \rightarrow North');
title('CS003, 10s integ, 40chans, 13:39');

%subplot(2, 2, 2);
figure(10);
imagesc(l, m, mean(sky_CS003_ts11_avg, 3).' .* mask);
set(gca, 'FontSize', 16, 'YDir', 'normal', 'XDir', 'reverse');
set(colorbar, 'FontSize', 16);
xlabel('East \leftarrow l \rightarrow West');
ylabel('South \leftarrow m \rightarrow North');
title('CS003, 10s integ, 40 chans, 13:49');

%subplot(2, 2, 3);
figure (11);
imagesc(l, m, mean(sky_CS003_ts12_avg, 3).' .* mask);
set(gca, 'FontSize', 16, 'YDir', 'normal', 'XDir', 'reverse');
set(colorbar, 'FontSize', 16);
xlabel('East \leftarrow l \rightarrow West');
ylabel('South \leftarrow m \rightarrow North');
title('CS003, 10s integ, 40chans, 13:59');

%subplot(2, 2, 4);
figure (12);
imagesc(l, m, mean(sky_CS003_ts13_avg, 3).' .* mask);
set(gca, 'FontSize', 16, 'YDir', 'normal', 'XDir', 'reverse');
set(colorbar, 'FontSize', 16);
xlabel('East \leftarrow l \rightarrow West');
ylabel('South \leftarrow m \rightarrow North');
title('CS003, 10s integ, 40 chans, 14:09');

%figure (2);
figure (13);
set(gcf, 'Position', [0, 500, 1120 840]);
%subplot(2, 2, 1);
imagesc(l, m, mean(sky_CS003_ts14_avg, 3).' .* mask);
set(gca, 'FontSize', 16, 'YDir', 'normal', 'XDir', 'reverse');
set(colorbar, 'FontSize', 16);
xlabel('East \leftarrow l \rightarrow West');
ylabel('South \leftarrow m \rightarrow North');
title('CS003, 10s integ, 40 chans, 14:19');

%subplot(2, 2, 2);
figure (14);
imagesc(l, m, mean(sky_CS003_ts15_avg, 3).' .* mask);
set(gca, 'FontSize', 16, 'YDir', 'normal', 'XDir', 'reverse');
set(colorbar, 'FontSize', 16);
xlabel('East \leftarrow l \rightarrow West');
ylabel('South \leftarrow m \rightarrow North');
title('CS003, 10s integ, 40 chans, 14:29');

%subplot(2, 2, 3);
figure (15);
imagesc(l, m, mean(sky_CS003_ts16_avg, 3).' .* mask);
set(gca, 'FontSize', 16, 'YDir', 'normal', 'XDir', 'reverse');
set(colorbar, 'FontSize', 16);
xlabel('East \leftarrow l \rightarrow West');
ylabel('South \leftarrow m \rightarrow North');
title('CS003, 10s integ, 40chans, 14:39');

%subplot(2, 2, 4);
figure (16);
imagesc(l, m, mean(sky_CS003_ts17_avg, 3).' .* mask);
set(gca, 'FontSize', 16, 'YDir', 'normal', 'XDir', 'reverse');
set(colorbar, 'FontSize', 16);
xlabel('East \leftarrow l \rightarrow West');
ylabel('South \leftarrow m \rightarrow North');
title('CS003, 10s integ, 40 chans, 14:49');

figure (17);
imagesc(l, m, mean(sky_CS003_ts18_avg, 3).' .* mask);
set(gca, 'FontSize', 16, 'YDir', 'normal', 'XDir', 'reverse');
set(colorbar, 'FontSize', 16);
xlabel('East \leftarrow l \rightarrow West');
ylabel('South \leftarrow m \rightarrow North');
title('CS003, 10s integ, 40 chans, 14:59');
