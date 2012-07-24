%function [cal,sigmas, Sigman]=track_statcal (acc, t_obs, freq, pos, srcsel, normal, restriction, maxrestriction, uvflag)
% Function prototype returning all iterations of the parameter estimation 
function [ghat, sigmahat, iter, cal, sigmas, Sigman, Sigmanhat]=track_statcal (acc, t_obs, freq, pos, srcsel, normal, restriction, maxrestriction, uvflag, maxiter, srclist3CR, sigmas0)

% [cal, sigmas, Sigman] =
%     statcal(acc, t_obs, freq, pos, srcsel, normal, restriction,
%             maxresitriction, uvflag)
%
% Wrapper function to provide the appropriate data to cal_ext
%
% arguments
% acc            : Nelem x Nelem x Nch array covariance matrix
% t_obs          : Nch x 1 vector of observing times, as JulianDay
% freq           : Nch x 1 vector of observing frequencies
% pos            : Nelem x 3 matrix with ITRF positions of the antennas
% srcsel         : Nsrc x 1 vector with source indices, use 0 for the Sun
% normal         : 3 x 1 normal vector to the station field
% restriction    : relative baseline restriction in wavelength
% maxrestriction : maximum absolute baseline restriction in meters
% uvflag         : Nelem x Nelem matrix to flag (reject) specific baselines
%
% return values
% cal    : Nelem x 1 vector with complex valued calibration corrections
% sigmas : Nsrc x 1 vector with estimated apparent source powers
% Sigman : Nelem x Nelem estimated array covariance matrix
%
% SJW, 2009
% modified on 18 May 2011 by SJW to use ITRF coordinates

% parameter section
diffstop = 1e-3;
c = 2.99792e8;
Nelem = size(pos, 1);
Nsb = length(freq);
% load srclist3CR % Now passing as a parameter to prevent reloads
nsrc = length(srcsel);
persistent first;
persistent flux;

% initialization
cal = zeros(Nsb, Nelem);
sigmas = zeros(Nsb, nsrc);
Sigman = zeros(Nsb, Nelem, Nelem);
u = meshgrid(pos(:, 1)) - meshgrid(pos(:, 1)).';
v = meshgrid(pos(:, 2)) - meshgrid(pos(:, 2)).';
w = meshgrid(pos(:, 3)) - meshgrid(pos(:, 3)).';
uvw = [u(:), v(:), w(:)];
uvdist = sqrt(sum(uvw.^2, 2) - (uvw * normal).^2);

for idx = 1:length(freq)
    % disp(['working on subband ' num2str(idx) ' of ' num2str(length(freq))]);
    % pause(0.01);
    % tic
    if (sum(srcsel == 0) ~= 0)
      %  [raSun, decSun] = SunRaDec(JulianDay(t_obs(idx)));          % J2000
	  %  disp (['JD ' num2str(JulianDay(t_obs(idx))) 'Sun RA/Dec' num2str(raSun) ' '  num2str(decSun)]);
       [raSun, decSun] = SunRaDec(t_obs);          % J2000
	   disp (['JD ' num2str(t_obs(idx)) 'Sun RA/Dec' num2str(raSun) ' '  num2str(decSun)]);
        rasrc = [srclist3CR(srcsel(srcsel ~= 0)).alpha, raSun].';   % B1950
        decsrc = [srclist3CR(srcsel(srcsel ~= 0)).delta, decSun].'; % B1950
        epoch = [true(length(srcsel(srcsel ~= 0)), 1); false];
    else
        rasrc = srclist3CR(srcsel).alpha;  % B1950
        decsrc = srclist3CR(srcsel).delta; % B1950
        epoch = true(length(srcsel), 1);
    end
   % srcpos = radectoITRF(rasrc, decsrc, epoch, JulianDay(t_obs(idx)))
    srcpos = radectoITRF(rasrc, decsrc, epoch, t_obs);
    % disp (['Source pos', num2str(srcpos)]);
    up = srcpos * normal > 0;
    disp ('statcal: Source coord. of visible srcs: '); disp (srcpos(up, :));
    A = exp(-(2 * pi * 1i * freq(idx) / c) * (pos * srcpos(up, :).'));
    Rhat = squeeze(acc(:, :, idx));

    mask = reshape(uvdist, [Nelem, Nelem]) < min([restriction * (c / freq(idx)), maxrestriction]);
    mask = mask | uvflag;

    if isempty(first)
        disp ('First time invocation: estimating model source flux using LSimaging..');
        first(1) = 0;
        flux = real(((abs(A' * A).^2) \ khatrirao(conj(A), A)') * (Rhat(:) .* (1 - mask(:))));
        flux = flux / flux(1);
        flux(flux < 0) = 0;
        disp ('Initial model flux: '); disp (flux'); % Size of flux'); size(flux)
        maxiter = 5;
    end 
    %[ghat, sigmahat, Sigmanhat] = cal_ext(Rhat, A, flux, mask, diffstop, maxiter);
    %iter = 1;
    %cal(idx, :) = conj(1./ghat);
    %sigmas(idx, up) = sigmahat ; %(:, iter);
    %Sigman(idx, :, :) = Sigmanhat; 

    % implementation of WALS method
    % estimate direction independent gains, apparent source powers and
    % receiver noise powers assuming that the source locations are correct
    [ghat, sigmahat, iter, Sigmanhat] = track_cal_ext(Rhat, A, flux, mask, diffstop, maxiter);
    flux = sigmahat (:, iter); % To pass on as initial estimate to the next dataset
    disp('Estimated flux for model sources: '); disp (flux); 
    
    cal(idx, :) = conj(1./ghat(:,iter));
    sigmas(idx, up) = sigmahat (:, iter);
    % Sigman(idx, :, :) = Sigmanhat;    
    %toc
end
