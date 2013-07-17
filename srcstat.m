function rhosrc = srcstat(srcflux, freq)

% Function to compute source density at the specified source flux at the
% specified frequency.
%
% The statistics are based on the analysis in Jaap Bregman's Ph.D. thesis
% for source statistics at frequencies below 1.4 GHz
%
% Arguments
% srcflux : vector with flux limits for which to compute the source density
%           of sources with at least that strength (in Jy)
% freq    : frequencies at which the limiting fluxes hold (in Hz)
%
% Return value
% rhosrc  : vector with source density for each limiting flux value
%
% SJW, 1 March 2012

% reference point data at 21 cm
Sref = [0.02e-3, 0.6e-3, 2e-3, 20e-3, 100e-3, 0.3, 1.7];       % Jy
rhoref = [4.12e7, 4.43e5, 1.88e5, 2.98e4, 5.92e3, 1.36e3, 81]; % sr^{-1}

% convert reference points to desired frequency
Sreff = zeros(length(Sref), 1);
if freq > 150e6
    Sreff(1)     = (1.4e9 / freq)^0.6 * Sref(1);
    Sreff(2)     = (1.4e9 / freq)^0.5 * Sref(2);
    Sreff(3)     = (1.4e9 / freq)^0.6 * Sref(3);
    Sreff(4:end) = (1.4e9 / freq)^0.8 * Sref(4:end);
else
    Sreff(1)     = (1.4e9 / 150e6)^0.6 * (150e6 / freq).^0.7 * Sref(1);
    Sreff(2)     = (1.4e9 / 150e6)^0.5 * (150e6 / freq).^0.7 * Sref(2);
    Sreff(3)     = (1.4e9 / 150e6)^0.6 * (150e6 / freq).^0.7 * Sref(3);
    Sreff(4)     = (1.4e9 / 150e6)^0.8 * (150e6 / freq).^0.7 * Sref(4);
    Sreff(5:end) = (1.4e9 / freq)^0.8 * Sref(5:end);
end

rhosrc = zeros(length(srcflux), 1);
rhosrc(srcflux < Sreff(1)) = NaN;
for idx = 2:length(Sreff)
    polcoeff = polyfit(log([Sreff(idx-1), Sreff(idx)]), log([rhoref(idx-1), rhoref(idx)]), 1);
    % modified to allow extrapolation to low and high flux values
    if idx == 2
        sel = srcflux < Sreff(idx);
    elseif idx == length(Sreff)
        sel = srcflux >= Sreff(idx-1);
    else
        sel = (srcflux >= Sreff(idx-1)) & (srcflux < Sreff(idx));
    end
    rhosrc(sel) = rhoref(idx-1) .* (srcflux(sel) / Sreff(idx-1)).^polcoeff(1);
end
