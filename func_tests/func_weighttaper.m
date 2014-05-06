% Script to explore the weighting and tapering parameter space for optimality.
% pep/04Apr14

%% LBA_OUTER
load ('../poslocal_outer.mat', 'poslocal', 'posITRF');
freq = 60000000;
flagant = [];
deb=1;



% Parameter space:
% Uniform weighting: Increasing cell radius
% Taper:
%   - Inc. inner taper, no outer
%   - Inc. outer taper, no inner
%   - Inc. inner and outer taper.

% Generate PSF without any taper/weight
fname = sprintf ('lba_outer_uni_notap_noweight');
tparm = []; wparm = []; gparm = [];
[l,m,psf,weight, intap, outtap] = genarraypsf ('poslocal.mat', flagant, freq, tparm, wparm, gparm, deb);
% Get figure handles
hand = get (0, 'Children');
print (hand(1), strcat (fname,'_psf.png'), '-dpng');
print (hand(2), strcat (fname,'_weights.png'), '-dpng');

fprintf (2, 'Increasing cell radius in uniform weighting, no tapering');
for ind=1:5:50
	% Generate weight parameters
	wparm.type = 'uniform'; wparm.cellrad = ind;
	fname = sprintf ('lba_outer_uni_notap_cr%02d', ind);
	[l,m,psf,weight, intap, outtap] = genarraypsf ('poslocal.mat', flagant, freq, [], wparm, [], [], [], deb);
	% Get figure handles
	hand = get (0, 'Children');
	print (hand(1), strcat (fname,'_psf.png'), '-dpng');
	print (hand(2), strcat (fname,'_weights.png'), '-dpng');
end;

fprintf (2, 'Increasing inner taper, no weighting, no outer taper.');
for ind=1:5:50
	% Generate taper parameters
	tparm.type = 'Gaussian';
	tparm.minlambda = 10;    % Spatial freq. cutoff lower limit.NOTE: Units of lambda.
	tparm.maxmeters = 350;	 % Spatial freq. cutoff upper limit. NOTE:Units of meters.
	tparm.pa(1) = ind; %0.2;	 % Inner taper gaussian sigx. NOTE: Units of lambda. 
	tparm.pa(2) = tparm.pa(1);% Inner taper gaussian sigy
	tparm.pa(3) = -1; 	     % Outer taper Gaussian sigx. NOTE: Units of lambda.
	tparm.pa(4) = tparm.pa(3); % Outer taper Gaussian sigy.

	% Generate weight parameters
	% wparm.type = 'uniform'; wparm.cellrad = ind;
	fname = sprintf ('lba_outer_intap_%02d_noouttap', ind);
	[l,m,psf,weight, intap, outtap] = genarraypsf ('poslocal.mat', flagant, freq, tparm, [], [], [], [], deb);
	% Get figure handles
	hand = get (0, 'Children');
	print (hand(1), strcat (fname,'_psf.png'), '-dpng');
	print (hand(2), strcat (fname,'_weights.png'), '-dpng');
end;

fprintf (2, 'Increasing outer taper, no weighting, no outer taper.');
for ind=1:10:350
	% Generate taper parameters
	tparm.type = 'Gaussian';
	tparm.minlambda = 10;    % Spatial freq. cutoff lower limit.NOTE: Units of lambda.
	tparm.maxmeters = 350;	 % Spatial freq. cutoff upper limit. NOTE:Units of meters.
	tparm.pa(1) = -1; %0.2;	 % Inner taper gaussian sigx. NOTE: Units of lambda. 
	tparm.pa(2) = tparm.pa(1);% Inner taper gaussian sigy
	tparm.pa(3) = ind; 	     % Outer taper Gaussian sigx. NOTE: Units of lambda.
	tparm.pa(4) = tparm.pa(3); % Outer taper Gaussian sigy.

	% Generate weight parameters
	% wparm.type = 'uniform'; wparm.cellrad = ind;
	fname = sprintf ('lba_outer_nowt_nointap_outtap_%02d', ind);
	[l,m,psf,weight, intap, outtap] = genarraypsf ('poslocal.mat', flagant, freq, tparm, [], [], [], [], deb);
	% Get figure handles
	hand = get (0, 'Children');
	print (hand(1), strcat (fname,'_psf.png'), '-dpng');
	print (hand(2), strcat (fname,'_weights.png'), '-dpng');
end;
