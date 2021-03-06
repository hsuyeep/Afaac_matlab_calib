% Program to test the iterative estimation of direction independent gains,
% after an initial solution has been estimated using a full fledged
% calibration. Operates on a single timeslice of correlated data
% pep/22Feb12
% Updated to operate on multiple timeslices, but just for the estimation of
% the model parameters.
% pep/15Mar12

%%
freq = 59951782.2;
C = 299792458; % speed of light
l = [-1:0.005:1]; m = l;
restriction = 4;             % Avoid uv values with <10 wavelengths length
maxrestriction = 60;          % Avoid uv values with >60 wavelengths length
fnames = {'SB000_1415-1420.bin', '1222_1227.bin'};
nfiles = 2;
Nelem = 288; 
nblines = Nelem * (Nelem + 1)/2; 
% FFT Imaging related
duv = 600/511;                % In meters, grid spacing in fft imager
Nuv = 512;                    % size of gridded visibility matrix
uvpad = 512;                  % specifies if any padding needs to be added
calmap = zeros (Nuv, Nuv);
calvis = calmap;

%%%%%%%%%%%%%%%%%% Variables for the unrolled versions of statcal and cal_ext%%%%%%%%%%%%%%%%%%
%% parameter section
diffstop = 1e-3;
maxiter = 10;
c = 2.99792e8;

%load srclist3CR
% initialization
load srclist3CR
srcsel =  [324, 283, 88, 179, 0];
nsrc = length(srcsel);
Nsb = length(freq);
cal = zeros(Nsb, Nelem);
sigmas = zeros(Nsb, nsrc);
Sigman = zeros(Nsb, Nelem, Nelem);
% Local horizon based coordinates of array in ITRF
load ('poslocal.mat', 'posITRF', 'poslocal'); 
normal = [0.598753, 0.072099, 0.797682].';
Nelem = size(posITRF, 1);

u = meshgrid(posITRF(:, 1)) - meshgrid(posITRF(:, 1)).';
v = meshgrid(posITRF(:, 2)) - meshgrid(posITRF(:, 2)).';
w = meshgrid(posITRF(:, 3)) - meshgrid(posITRF(:, 3)).';
uvw = [u(:), v(:), w(:)];
uvdist = sqrt(sum(uvw.^2, 2) - (uvw * normal).^2);
    
% Generate uv coordinates in local horizon coordinate system, needed for imaging
uloc = meshgrid (poslocal(:,1)) - meshgrid (poslocal (:,1)).';
vloc = meshgrid (poslocal(:,2)) - meshgrid (poslocal (:,2)).';

%%

fid = fopen (fnames{1}, 'rb');
t_obs1 = fread (fid, 1, 'double');
disp (['Timeslice: ' num2str(t_obs1, '%f') ' (MJD secs) for initial estimation']);
% Convert back to Julian day, as desired by track_statcal. NOTE:
% JulianDay () should not be called now!
t_obs = t_obs1/86400 + 2400000.5;

% Reading real and imaginary, available as a stream of floats.
% even floats being real parts, odd floats being imag
a = fread (fid, 2*nblines, 'float'); % Read one ccm worth
comp = complex (a(1:2:length(a)), a(2:2:length(a))); % to complex
% create instantaneous ccm from vector
acm = triu (ones (Nelem));
acm (acm == 1) = comp;
acc = acm + acm' - diag(diag(acm));
uvflag = eye (Nelem);
 
% Carry out a first level estimation of all parameters, using a full
% fledged calibration on the first, single timeslice, in order to generate
% initial estimates for direction independent gains, apparent source
% powers, and system noise.
[cal0, sigmas0, Sigman0] = statcal(acc, t_obs, freq, posITRF, srcsel, normal, 4, 20, eye(Nelem));
ghat0 = conj (1./cal0);

disp ('Initial estimates of model src flux:'); disp (sigmas0); 
disp (['Size of sigmas: ' num2str(length(sigmas0)) ', min/max noise: ' num2str(min(min(abs(squeeze (Sigman0))))) '/' num2str(max(max(abs(squeeze (Sigman0)))))]);

prevgainwin = figure ();
subplot (2,1,1);
plot (abs(ghat0)); title ('Prev iteration amplitude solution');
subplot (2,1,2);
plot (angle (ghat0)); title ('Prev iteration phase solution');

currgainwin = figure ();
% Iterate over several timeslices, using each for a single iteration while
% estimating model parameters.

%% while ~feof (fid)
 for i = 1:10 % Initially only for 10 timeslices.
    % Read in another timeslice to work with
    t_obs1 = fread (fid, 1, 'double');
    disp (['Timeslice: ' num2str(t_obs1, '%f') ' (MJD secs)']);
    % Convert back to Julian day, as desired by track_statcal. NOTE:
    % JulianDay () should not be called now!
    t_obs = t_obs1/86400 + 2400000.5;
    a = fread (fid, 2*nblines, 'float'); % Read one ccm worth
    comp = complex (a(1:2:length(a)), a(2:2:length(a))); % to complex
    % create instantaneous ccm from vector
    acm = triu (ones (Nelem));
    acm (acm == 1) = comp;
    acc = acm + acm' - diag(diag(acm));
    idx = 1;

    %%
    for idx = 1:length(freq)
        % disp(['working on subband ' num2str(idx) ' of ' num2str(length(freq))]);
        % tic
        if (sum(srcsel == 0) ~= 0)
            % [raSun, decSun] = SunRaDec(JulianDay(t_obs(idx)));          % J2000
            [raSun, decSun] = SunRaDec(t_obs(idx));          % J2000
            rasrc = [srclist3CR(srcsel(srcsel ~= 0)).alpha, raSun].';   % B1950
            decsrc = [srclist3CR(srcsel(srcsel ~= 0)).delta, decSun].'; % B1950
            epoch = [true(length(srcsel(srcsel ~= 0)), 1); false];
        else
            rasrc = srclist3CR(srcsel).alpha;  % B1950
            decsrc = srclist3CR(srcsel).delta; % B1950
            epoch = true(length(srcsel), 1);
        end
        % srcpos = radectoITRF(rasrc, decsrc, epoch, JulianDay(t_obs(idx)));
        srcpos = radectoITRF(rasrc, decsrc, epoch, t_obs);

        up = srcpos * normal > 0;
        disp (['Number of model sources above the horizon: ' num2str(sum(up))]);
        A = exp(-(2 * pi * 1i * freq(idx) / c) * (posITRF * srcpos(up, :).'));
        Rhat = squeeze(acc(:, :, idx)); % Observed acm

        mask = reshape(uvdist, [Nelem, Nelem]) < min([restriction * (c / freq(idx)), maxrestriction]);
        mask = mask | uvflag;

        % NOTE: Changed the 'flux' variable name to sigmas, to match cal_ext ()
        % usage. Should this be done? as it is already done in the first
        % initialization call to statcal(). Maybe we should just use
        % sigmas0...
        sigmas = real(((abs(A' * A).^2) \ khatrirao(conj(A), A)') * (Rhat(:) .* (1 - mask(:))));
        sigmas = sigmas / sigmas(1);
        sigmas(sigmas < 0) = 0;
        disp ('Statcal: model src flux estimates from LSimaging: ');
        disp (sigmas);


    %% implementation of WALS method
    % estimate direction independent gains, apparent source powers and
    % receiver noise powers assuming that the source locations are correct
    % Function call as seen in statcal:
    %    [ghat, sigmahat, Sigmanhat] = cal_ext(Rhat, A, flux, mask, diffstop,
    %    maxiter);
    % Function parameter from cal_ext:
    %   [g, sigmas, Sigma_n] = cal_ext (Rhat, A, sigmas, mask, varargin);

    % cal_ext code, after unwrapping from function
    % parameters
    [Nelem, Nsrc] = size(A);


    % initialization
    % ghat = zeros(Nelem, maxiter+1);
    ghat1 = zeros (Nelem, 1);
    sigmahat = zeros(Nsrc, maxiter+1);
    sigmas(~isfinite(sigmas)) = 1; % NOTE: Can use estimate of sigma from data, and not from LS imaging?
    sigmahat(:, 1) = sigmas;
    % Should this be done as well? Can use the noise estimate from the
    % previous iteration itself...
    % Sigma_n = zeros(Nelem);

    % Objects to be created for carrying out the iterative, direction
    % independent gain estimation procedure, with the iterations being
    % different timeslices.
    sigmas_mat = diag (sigmas);     % Include dir. dependent gains in dir. of every model source
    G = diag (conj (1./cal0));      % dir. independent gain est. from earlier timeslice.
    R0 = A * sigmas_mat * A' .* (1 - mask); % Model array response matrix
    Rhat_mask = Rhat .* (1 - mask); % Observed data, after removal of unwanted visibilities
    Rhat_clean = Rhat_mask - squeeze (Sigman0); % After removing diagonal noise terms

    %% implementation using WALS, using only a single iteration on the
    % current timeslice, but with initial solutions available from
    % previous timeslices.
    iter = 1;
    % for iter = 2:maxiter+1
        %%
        % estimate g using baseline restriction
        % tic
        % Implementation of the iterative direction independent gain estimate,
        % as per the 'Multisource selfcal paper'. NOTE: noise estimate is also
        % taken from the previous timeslice. Also, the paper mentions that this
        % method may not converge to the global minimum if the initial estimate
        % is bad, but our estimate is via a full fledged iterative solution for
        % g_est.
        R = G * R0 * G' + squeeze (Sigman0);
        rinv = inv (R);
        rinv_G_r0 = rinv * G * R0;
        % g_est = inv (conj (R0' * G' * rinv_G_r0) .* rinv) * khatrirao (conj (rinv_G_r0), rinv)' * Rhat_clean(:); 
        ghat1 = inv (conj (R0' * G' * rinv_G_r0) .* rinv) * khatrirao (conj (rinv_G_r0), rinv)' * Rhat_clean(:); 

        %rinv =  R \ eye(Nelem);  % NOTE: Was doing this, as it is expected to be faster, but was actually slower 
        %t1 = (conj (R0' * G' * rinv * G * R0) .* rinv) \  khatrirao (conj (R \ G * R0), rinv)' * Rhat_clean (:);
        %t2 = khatrirao (conj (R \ G * R0), rinv)';
        %g_est = t1\t2 * Rhat_clean (:);
        %g_est = t1 * Rhat_clean (:);
        % toc;

    %     alpha = computeAlphaW(Rhat .* (1 - mask), (A * diag(sigmahat(:, iter-1)) * A') .* (1 - mask));
    %     toc
    %     tic
    %     if (sum(~isfinite(alpha(:))) == 0)
    %         [v, d] = eig(alpha);
    %         [~, idx] = max(diag(d));
    %         ghat(:, iter) = conj(v(:, idx));
    %     else
    %         ghat(:, iter) = 1;
    %     end
    %     
    %     GA = diag(ghat(:, iter)) * A;
    %     Rest = GA * diag(sigmahat(:, iter-1)) * GA';
    %     compidx = (1 - mask) ~= 0;
    %     normg = abs(sqrt(pinv(Rest(compidx)) * Rhat(compidx)));
    %     ghat(:, iter) = normg * ghat(:, iter) / (ghat(1, iter) / abs(ghat(1, iter)));
    %     ghat(~isfinite(ghat(:, iter)), iter) = 1;
    %    toc

        %% estimate sigmahat using sigmanhat and ghat
        % tic
        invR = inv(Rhat);
        % GA = diag(ghat(:, iter)) * A; % new line, use normalized G
        % GA = diag (g_est) * A;
        GA = diag (ghat1) * A;
        % sigmahat(:, iter) = real(inv(abs(conj(GA' * invR * GA)).^2) *
        % diag(GA' * invR * (Rhat - Sigma_n) * invR * GA));
        % Now using the initial estimate of system noise.
        sigmahat(:, iter) = real(inv(abs(conj(GA' * invR * GA)).^2) * diag(GA' * invR * (Rhat - squeeze(Sigman0)) * invR * GA));
        if sum(~isfinite(sigmahat(:, iter))) ~= 0
            sigmahat(:, iter) = sigmahat(:, iter-1);
        end
        % use first source as amplitude reference
        if (sigmahat(1, iter) ~= 0)
            sigmahat(:, iter) = sigmahat(:, iter) / sigmahat(1, iter);
        end
        % remove negative values
        sigmahat(:, iter) = max(sigmahat(:, iter), 0);
        % toc

        % estimate new sigmanhat using sigmahat and ghat.
        % NOTE: To plot this (e.g. to cmp against initial noise estimate),
        % imagesc (abs (squeeze (Sigma_n)) - abs (diag(diag(squeeze (Sigma_n)))));
        Sigman1 = (Rhat - GA * diag(sigmahat(:, iter)) * GA') .* mask;

        % test for convergence
        % theta_prev = [ghat(:, iter-1); sigmahat(:, iter-1)];
        % theta = [ghat(:, iter); sigmahat(:, iter)];

        % theta_prev = [g_est; sigmahat(:, iter-1)];
        % theta = [g_est; sigmahat(:, iter)];
        % cost = abs(pinv(theta_prev) * theta - 1);
        % if (abs(pinv(theta_prev) * theta - 1) < diffstop)
        % if (cost < diffstop)
        %    break
        % end

        % disp (['iter: ' num2str(iter) 'cost :' num2str(cost)]);
        
    %end
    figure (currgainwin);    
    subplot (2,1,1);
    plot (abs(ghat1));
    title (['Current iteration amp. solutions. Time: ' num2str(t_obs, '%f')]);
    subplot (2,1,2);
    plot (angle(ghat1));
    title (['Current iteration phase solutions. Time: ' num2str(t_obs, '%f')]);
    
    
    disp(['optimization stopped after ' num2str(iter) ' iteration(s).']);
    % g = ghat(:, iter);
    sigmas = sigmahat(:, iter);

    
    
    % From statcal
    % cal(idx, :) = conj(1./ghat (:, iter)); % NOTE: Needed to change this for picking up last iteration
    %sigmas(idx, up) = sigmahat; % Returning sigmahat
    %Sigman(idx, :, :) = Sigmanhat;    
    %toc
    end
 end


