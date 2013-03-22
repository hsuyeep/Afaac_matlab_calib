% fitAmplModel
%
% This script is used by the script makeCalTable.m to fit a gain amplitude
% model to the station calibration solutions. Currently, the model is simply a
% constant amplitude over frequency. This script includes outlier removal and
% uses the median function for additional robustness. It also displays a
% messages stating how many points are used. A low number of data points (a few
% hundred or less) generally indicates a problem such as signal path failures.
%
% SJW, December 2010
% modified by SJW, October 2011

% loop over all unprocessed signal paths
abscal = zeros(Nant, 1);
Nsample = zeros(Nant, 1);
for antidx = 1:Nant
    % collect data for this antenna and remove flagged data
    abscalant = squeeze(abs(cal(:, :, antidx)));
    data = abscalant(abscalant(:) ~= 0 & abscalant(:) ~= 1);
    % compute average gain using median
    abscal(antidx) = median(data);
    % use this first estimate to remove outliers
    meanabsdiff = mean(abs(data - abscal(antidx)));
    data = data(abs(data - abscal(antidx)) < 3 * meanabsdiff);
    Nsample(antidx) = length(data);
    % determine average gain using median function
    abscal(antidx) = median(data);
end

disp(['The amplitude fits are on average based on ' num2str(mean(Nsample)) ' data points']);
disp(['with standard deviation ' num2str(std(Nsample)) '.']);
badant = 1:Nant;
badant = badant(Nsample < 0.2 * mean(Nsample));
if (~isempty(badant))
    disp('The following antennas have fewer than 20% of the average data points:');
    disp(num2str(badant));
end
