% Script to generate single timeslice calibration efficiency.
% pep/31Jan14

% Open the uncalibrated visibilities
fid = fopen ('~/WORK/AARTFAAC/Reobs/11Jul12/LBA_OUTER_BAND_SPREAD/SB002_LBA_OUTER_SPREAD_1ch.bin', 'rb');
flagant = [51, 239, 273];

load ('poslocal_outer.mat', 'poslocal', 'posITRF'); % Load baseline coordinates

% Read in a timeslice
[acc, tacc, freq] = readms2float (fid, -1, -1, 288); 
fprintf (2, 'Time: %f, Freq: %f\n', tacc, freq);
sol = pelican_sunAteamsub (acc, tacc, freq, eye(288), flagant, 0, 1, [], [], 'poslocal_outer.mat');

% Handle flagging.
normal  = [0.598753, 0.072099, 0.797682].'; % Normal to CS002
antmask = zeros (size (acc));
posmask = zeros (size (posITRF));
rem_ants = length(acc) - length(flagant);

for ind = 1:length(flagant)
	antmask (flagant(ind), :) = 1; antmask (:,flagant(ind)) = 1;
	posmask (flagant(ind), :) = 1;
end
posITRF_fl = reshape(posITRF(posmask ~=1), [rem_ants, 3]);
u=meshgrid(posITRF_fl(:, 1))-meshgrid(posITRF_fl(:, 1)).';
v=meshgrid(posITRF_fl(:, 2))-meshgrid(posITRF_fl(:, 2)).';
w=meshgrid(posITRF_fl(:, 3))-meshgrid(posITRF_fl(:, 3)).';
uvw = [u(:), v(:), w(:)];
uvdist = sqrt(sum(uvw.^2, 2) - (uvw * normal).^2);

% Generate visibilities from model sky for this timeslice (Generated from within 
% pelican_sunAteamsub.m for this timeslice.
load ('4848774150_Ateammodelvis.mat');

% Generate visibility ratio with model 
visrat = sol.calvis(:) ./ RAteam(:);

% Arrange visibilities in increasing order of visibility power
visratpwr = abs(visrat);
[sortvisratpwr, sortind] = sort (visratpwr, 'ascend');

 % plot residual visibility phase as a function of increasing visibility power.
plot (

