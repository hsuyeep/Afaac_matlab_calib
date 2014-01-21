% Script to test out the visibility tapering function.
% pep/20Aug13
function func_tapertest (minlambda, maxmeters, insig, outsig)

%%%%%%  Part dealing with calling the taper function %%%%%%%%
		posfilename = 'poslocal_outer.mat';
		load (posfilename, 'posITRF', 'poslocal');
		Nelem = length (posITRF);

		parm.type = 'Gaussian';
	    parm.minlambda = minlambda; % lambda 
		parm.maxmeters = maxmeters; % meters 
		parm.pa(1) = insig;
		parm.pa(2) = insig;
		parm.pa(3) = outsig;
		parm.pa(4) = outsig;

		[intap, outtap, den, mask, uvdist] = taper (posITRF, parm, -1, 60e6, 4);
%		figure;
%		subplot (121); imagesc (mask); 
%		subplot (122); plot (uvdist(:), (mask(:)), '.');

%%%%%%%% Part dealing with analysis of calibration solutions as a function of
%%%%%%%% spatial filter extent. Assumes the presence of the variable 'sols', whic
%% is an array of structures containing the solutions as returned by 
%% pelican_sunAteamsub.m
%
%% show all gain solutions
%for ind = 1:16
%	plot (abs(sols(ind).gainsol));
%	title (sprintf ('%d m cutoff', ind*10));
%	pause (0.5);
%end
%
%for ind = 1:16
%	plot (angle(sols(ind).gainsol));
%	title (sprintf ('%d m cutoff', ind*10));
%	pause (0.5);
%end;
%
%for ind = 1:16
%	plot (abs(diag(sols(ind).sigman)));
%	title (sprintf ('%d m cutoff', ind*10));
%	pause(0.5);
%end
%
%% Errors relative to 50m taper
%for ind = 1:16
%	plot (angle(sols(ind).gainsol) - angle(sols(5).gainsol), '.');
%	title (sprintf ('%d m cutoff', ind*10));
%	pause
%end
%
%for ind = 1:16
%	plot (abs(sols(ind).gainsol - sols(5).gainsol));
%	title (sprintf ('%d m cutoff', ind*10));
%	pause;
%end;
%
%for ind = 1:16
%	plot (abs(sols(ind).sigman - sols(5).sigman));
%	title (sprintf ('%d m cutoff', ind*10));
%	pause;
%end
%
%% plot number of iterations
%plot ([10:10:160], [sols(:).stefcal_iters],'d')
%xlabel ('Spatial cutoff (m)') ylabel ('Stefcal iterations');
%title ('Variation of stefcal iterations with spatial filter cutoff');
%
%plot ([10:10:160], [sols(:).pinv_sol],'d');
%xlabel ('Spatial cutoff (m)') ylabel ('Stefcal iterations');
%title ('Variation of stefcal iterations with spatial filter cutoff');
%
%plot ([10:10:160], [sols(:).stefcal_tol],'d');
%xlabel ('Spatial cutoff (m)') ylabel ('Stefcal iterations');
%title ('Variation of stefcal iterations with spatial filter cutoff');
%
%% Plot relative sigman contribution
%for ind = 1:16
%	mask = (sols (ind).sigman == 0);
%	plot (ind, sum(sum(abs(sols(ind).sigman)))/sum(sum(abs(acc(mask == 1)))), '*');
%	hold on;
%end
%grid on
%xlabel ('Spatial cutoff (x10m)');
%ylabel ('Integrated sigman amp frac. of uncal flux');
%title ('Estimated integ. sigman flux as fraction of uncal flux');
%
%
%% sigmas variation
%for ind = 1:16
%	plot (ind, sols(ind).sigmas(2) - sols(5).sigmas(2), '*');
%	hold on;
%	pause;
%end;
%
%% Cumulative distribution of power as a function of baseline length.
[sort_uvdist, sort_ind] = sort(uvdist, 1, 'ascend');
b = abs(acc(:))/sum(abs(acc(:))); % Normalize visibility fluxes to 1, for %
plot (uvdist (sort_ind), cumsum(b(sort_ind)),'.')
