% Program to implement a 'tracking calibration', essentially by combining
% SJW's statcal.m and cal_ext.m. Further, only a single iteration is
% performed in cal_ext.m, but with the facility of passing an initial
% estimate of estimated complex gain, complex noise and estimated source
% flux. 
% NOTE: Assumess that visibility mask and uvw is provided.
function [cal1, sigmas1, Sigman1] = tracking_cal (acc, t_obs, freq, pos, srcsel, srclist,  normal, mask) %, cal0, sigmas0, Sigman0, uvw)

c = 2.99792e8;
Nelem = size (pos, 1);
nsrc = length (srcsel);

% for cal_ext.m
ghat = zeros (Nelem, 1);

cal    = zeros(Nelem);
sigmas = zeros(nsrc);
Sigman = zeros(Nelem, Nelem);

if (sum(srcsel == 0) ~= 0)
    [raSun, decSun] = SunRaDec(JulianDay(t_obs));          % J2000
    rasrc  = [srclist(srcsel(srcsel ~= 0)).alpha, raSun].';   % B1950
    decsrc = [srclist(srcsel(srcsel ~= 0)).delta, decSun].'; % B1950
    epoch  = [true(length(srcsel(srcsel ~= 0)), 1); false];
else
    rasrc  = srclist(srcsel).alpha;  % B1950
    decsrc = srclist(srcsel).delta; % B1950
    epoch  = true(length(srcsel), 1);
end
    srcpos = radectoITRF(rasrc, decsrc, epoch, JulianDay(t_obs));
    up = srcpos * normal > 0;
    A = exp(-(2 * pi * 1i * freq / c) * (pos * srcpos(up, :).'));
    Rhat = squeeze(acc(:, :));

    % mask = reshape(uvdist, [Nelem, Nelem]) < min([restriction * (c / freq(idx)), maxrestriction]);
    % mask = mask | uvflag;

    flux = real(((abs(A' * A).^2) \ khatrirao(conj(A), A)') * (Rhat(:) .* (1 - mask(:))));
    flux = flux / flux(1);
    flux(flux < 0) = 0;

%     % implementation of WALS method, from cal_ext.m
%     %ghat = zeros(Nelem, maxiter+1);
%     sigmahat = zeros(Nsrc, maxiter+1);
%     sigmas(~isfinite(sigmas)) = 1;
%     sigmahat(:, 1) = sigmas;
%     Sigma_n = zeros(Nelem);
% 
% % implementation using WALS
% maxiter = 10;
% for iter = 2:maxiter+1
%     
%     % estimate g using baseline restriction
%     alpha = computeAlphaW(Rhat .* (1 - mask), (A * diag(sigmahat(:, iter-1)) * A') .* (1 - mask));
%     if (sum(~isfinite(alpha(:))) == 0)
%         [v, d] = eig(alpha);
%         [~, idx] = max(diag(d));
%         ghat(:, iter) = conj(v(:, idx));
%     else
%         ghat(:, iter) = 1;
%     end
%     GA = diag(ghat(:, iter)) * A;
%     Rest = GA * diag(sigmahat(:, iter-1)) * GA';
%     compidx = (1 - mask) ~= 0;
%     normg = abs(sqrt(pinv(Rest(compidx)) * Rhat(compidx)));
%     ghat(:, iter) = normg * ghat(:, iter) / (ghat(1, iter) / abs(ghat(1, iter)));
%     ghat(~isfinite(ghat(:, iter)), iter) = 1;
% 
%     % estimate sigmahat using sigmanhat and ghat
%     invR = inv(Rhat);
%     GA = diag(ghat(:, iter)) * A; % new line, use normalized G
%     sigmahat(:, iter) = real(inv(abs(conj(GA' * invR * GA)).^2) * diag(GA' * invR * (Rhat - Sigma_n) * invR * GA));
%     if sum(~isfinite(sigmahat(:, iter))) ~= 0
%         sigmahat(:, iter) = sigmahat(:, iter-1);
%     end
%     % use first source as amplitude reference
%     if (sigmahat(1, iter) ~= 0)
%         sigmahat(:, iter) = sigmahat(:, iter) / sigmahat(1, iter);
%     end
%     % remove negative values
%     sigmahat(:, iter) = max(sigmahat(:, iter), 0);
%     
%     % estimate sigmanhat using sigmahat and ghat
%     Sigma_n = (Rhat - GA * diag(sigmahat(:, iter)) * GA') .* mask;
% 
%     % test for convergence
%     theta_prev = [ghat(:, iter-1); sigmahat(:, iter-1)];
%     theta = [ghat(:, iter); sigmahat(:, iter)];
%     if (abs(pinv(theta_prev) * theta - 1) < diffstop)
%         break
%     end
% end
% disp(['optimization stopped after ' num2str(iter) ' iteration(s).']);
% g = ghat(:, iter);
% sigmas = sigmahat(:, iter);

    % estimate direction independent gains, apparent source powers and
    % receiver noise powers assuming that the source locations are correct
    diffstop = 1e-3;
    maxiter = 10;
    [ghat, sigmahat, Sigmanhat] = cal_ext(Rhat, A, flux, mask, diffstop, maxiter);
    cal1 = conj(1./ghat);
    sigmas1(up) = sigmahat;
    Sigman1 = Sigmanhat;    
    
    

