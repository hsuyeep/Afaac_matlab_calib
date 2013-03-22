% fitAmplModel
%
% This script is used by the script makeCalTable.m to fit a phase mode to the
% gain solutions. Currently, the phase model assumes a linear frequency
% response. This script first determines the average phase for each frequency
% over all available snaspshots in time using the median function for
% robustness to outliers. By setting the antstart and sbflag variables in the
% Matlab workspace, specific subbands can be flagged for a specific antenna in
% case phase unwrapping does not work correctly.
%
% SJW, Decemeber 2010
% modified by SJW, October 2011

% loop over all unprocessed antennas
argcaldata = angle(cal);
argcal = zeros(Nant, Nfreq);
plotdata{16} = [];
plotidx{16} = [];
figure
set(gcf, 'Position', [0, 0, 1280 1024]);
for antidx = 1:Nant
    % collect one data point for each subband by taking average over time
    data = zeros(1, Nfreq);
    for fidx = 1:Nfreq
        sbdata = argcaldata(:, fidx, antidx);
        sbdata = sbdata(sbdata ~=0);
        data(fidx) = mean(exp(1i * sbdata));
    end
    data(abs(data) <= 0.98) = NaN;
    data(~isnan(data)) = angle(data(~isnan(data)));

    % check for smoothness of the data (no steep gradients)
    testdata = data(isfinite(data));
    unsmooth = (testdata - [testdata(2:end), NaN] > pi/4) | (testdata - [NaN, testdata(1:end-1)] > pi/4);
    testdata(unsmooth) = NaN;
    data(isfinite(data)) = testdata;
    
    % apply user specified flags and unwrap the data
    data(sbflag{antidx}) = NaN;
    data = unwrap(data);
    sel = isfinite(data); 

    % try to fit a linear model
    coeff = polyfit(freq(sel), data(sel), 1);
    argcal(antidx, :) = polyval(coeff, freq);
    % remove outliers and make a new fit
    meanabsdiff = mean(abs(data(sel) - argcal(antidx, sel)));
    fsel = freq(sel);
    cleandata = data(sel);
    sel2 = abs(cleandata - argcal(antidx, sel)) < 3 * meanabsdiff;
    coeff = polyfit(fsel(sel2), cleandata(sel2), 1);
    argcal(antidx, :) = polyval(coeff, freq);
    
    % plot the result for check by the user
    if (mod(antidx, 16) ~= 0)
        plotdata{mod(antidx, 16)} = cleandata(sel2);
        sbidx = 1:512;
        sbidx = sbidx(sel);
        plotidx{mod(antidx, 16)} = sbidx(sel2);
    else
        sbidx = 1:512;
        sbidx = sbidx(sel);
        plotidx{16} = sbidx(sel2);
        plotdata{16} = cleandata(sel2);
        for idx = 1:16
            subplot(4, 4, idx)
            plot(1:512, argcal(antidx - 16 + idx, :), '-', plotidx{idx}, plotdata{idx}, '.');
            title(['phase fit for antenna ' num2str(antidx - 16 + idx)]);
        end
        pause;
    end
end
close