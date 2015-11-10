% Script to compute the flux of calibrator sources, based on Anna Scaife's models.
% pep/06Nov15
% Arguments:
%   srcname :   Name of calibrator whose flux is required.
%   freq    :   Frequency (Hz) at which the flux is required.

function flux = gencalflux (srcname, freq)
    % Fit parameters from Scaife's model fits
    fitparms = containers.Map ({'3C48', '3C147', '3C196', '3C286'}, { [ 64.768  -0.387  -0.420   0.181], [ 66.738  -0.022  -1.012   0.549], [ 83.084  -0.699  -0.110   0    ], [ 27.477  -0.158   0.032  -0.180] });

    assert (fitparms.isKey(upper (srcname)));

    % Model specification
    coe = log10(freq/150e6);
    mod = [1 coe  coe^2 coe^3];

    flux = fitparms(upper(srcname)).*mod;
