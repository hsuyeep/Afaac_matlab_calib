function AoverT = findAoverT(cal, Sigman, freq, t_obs)

% AoverT = findAoverT(cal, Sigman, freq, t_obs)
%
% Simple estimate of A/T towards Cas A using the method described in [1] and
% the Cas A flux from [2]. Note that, as indicated in [1], improved Cas A
% flux estimates are available in [3] and [4], so the value computed here
% should only be used as an initial assessment of station performance.
%
% Arguments
% cal    : Nsb x Nelem matrix with element based complex gain corrections
% Sigman : Nsb x Nelem x Nelem cube of estimated noise cov. matrices
% freq   : Nsb x 1 vector containing subband center frequencies
% t_obs  : Nsb x 1 vector containing data acquisition time for each subband
%
% Return value
% AoverT : Nsb x 1 vector containing A/T estimate per element for each subband
%
% SJW, December 2010
%
% [1] Stefan J. Wijnholds and Wim A. van Cappellen, "In Situ Antenna
%     Performance Estimation of the LOFAR Phased Array Radio Telescope", IEEE
%     Transactions on Antennas and Propagation, 2011, in press
% [2] J.W.M. Baars et al., "The Absolute Spectrum of Cas A: An Accurate Flux
%     Density Scale and a Set of Secondary Calibrators", Astronomy & Astrophy-
%     sics, V61, pp99-106, 1977
% [3] Daniel E. Reichart and Andrew W. Stephens, "The Fading of Supernova
%     Remnant Cassiopeia A from 38 MHz to 16.5 GHz from 1949 to 1999 with New
%     Observations at 1405 MHz", The Astrophysical Journal, no. 537, pp904-908,
%     2000
% [4] J.F. Helmboldt and N.E. Kassim, "The Evolution of Cas A at Low Radio
%     Frequencies", The Astronomical Journal, V138, pp838-844, July 2009

% Initialization
k = 1.38066e-23;
Nfreq = length(freq);
PCas = zeros(Nfreq, 1);
Pac = zeros(Nfreq, 1);

% Determine measured power on Cas A (assumed to be the flux reference) and
% the autocorrelation power for each frequency channel
for ch = 1:Nfreq
    PCas(ch) = mean(abs(1 ./ cal(ch, :)))^2;
    Pac(ch) = mean(real(diag(squeeze(Sigman(ch, :, :)))));
end

% Compute Cas A flux following [2]
SCas1965 = 10.^(5.625 - 0.634 * log10(freq / 1e6) - 0.023 * log10(freq / 1e6).^2);
decay = (0.97 - 0.30 * log10(freq / 1e9)) / 100;
SCas = 1e-26 * SCas1965 .* (1 - decay).^((t_obs - datenum(1965, 1, 1, 0, 0, 0)) / 365.24);

% Compute A/T
AoverT = (PCas ./ Pac) .* (2 * k ./ SCas);
