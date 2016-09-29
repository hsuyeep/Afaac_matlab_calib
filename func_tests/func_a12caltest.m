% Script to test the existing calibration mechanism on A-12 data.
% pep/08Sep16

load ('/dop312_0/prasad/GPU_CORR_DAT/afaac-12/lbaouter_06Nv15/sb0_1447252930.vis_1447252960-1447252969.mat');
fobs = 195312.5*296; % Hz.
addpath '~/Documents/AARTFAAC_Project_SW_system_plan/afaac_GPU_interface/src'

% Looked at the visi. ampl in log space, and determined the following flags:
flagant = [324, 373, 419, 492, 527, 537, 538];
[acm_t, tmjdsec,fobs,map,l] = gengpuimg (acm, 576,tobs,fobs,[1:15],[],[],[],0,0);

acm_uncal= zeros (576, 576);
acm_uncal(:, :) = conj (acm_t(1,:,:,1));

uvflag = eye(576);

% Imaging related setup
gridparm.type = 'pillbox';
gridparm.lim  = 0;

% Control imaged field of view, independently of frequency.
% obs.gridparm.duv  = (0.5/180)*obs.fov;
gridparm.duv  = 0.5;

gridparm.Nuv  = 2048;    % size of gridded visibility matrix
gridparm.uvpad= 2048;    % specifies if any padding needs to be added
gridparm.fft  = 1;
load ('poslocal_afaac12_outer.mat');

% Taper function parameters
tparm.type = 'Gaussian';
tparm.minlambda = 10;
tparm.maxmeters = 1200;   % NOTE: Units of meters.
tparm.pa(1) = 0;         % NOTE: Units of lambda.
tparm.pa(2) = tparm.pa(1); % Inner taper sigx/sigy
tparm.pa(3) = -1;    % NOTE: Units of lambda.
tparm.pa(4) = tparm.pa(3); % Outer taper sigx/sigy

% Weighting function parameters
wparm.type = 'Uniform';
wparm.cellrad = 10;

% [l, m, psf, weight, intap, outtap,uvdist, vispad] = genarraypsf ('poslocal_afaac12_outer.mat', flagant, fobs,  tparm, [], gridparm, 1);
% Fix the coordinates of the mixed up stations
% Emperically determined that:
% CS011 uses the coords of CS001 (put stn 09 coords in stn 07)
% CS013 uses the coords of CS011 (put stn 07 coords in stn 08)
% CS001 uses the coords of CS013 (put stn 08 coords in stn 09)
% NOW FIXED IN poslocal_afaac12_outer.mat
%stn07 = poslocal (6*48+1:7*48,:);
%stn08 = poslocal (7*48+1:8*48,:);
%stn09 = poslocal (8*48+1:9*48,:);
%poslocal (6*48+1:7*48,:) = stn09;
%poslocal (7*48+1:8*48,:) = stn07;
%poslocal (8*48+1:9*48,:) = stn08;
%
normal  = [0.598753, 0.072099, 0.797682].'; % Normal to CS002
uloc = meshgrid (poslocal(:,1)) - meshgrid (poslocal (:,1)).';
vloc = meshgrid (poslocal(:,2)) - meshgrid (poslocal (:,2)).';
wloc = meshgrid (poslocal(:,3)) - meshgrid (poslocal (:,3)).';
[uloc_a12, vloc_a12, wloc_a12] = gen_flagged_uvloc (uloc, vloc, wloc, flagant);
uvw_a12 = [uloc_a12(:), vloc_a12(:), wloc_a12(:)];
uvdist_a12 = sqrt(sum(uvw_a12.^2, 2) - (uvw_a12 * normal).^2);

[uloc_a6, vloc_a6, wloc_a6] = gen_flagged_uvloc (uloc, vloc, wloc, union([289:576], flagant));
uvw_a6 = [uloc_a6(:), vloc_a6(:), wloc_a6(:)];
uvdist_a6 = sqrt(sum(uvw_a6.^2, 2) - (uvw_a6 * normal).^2);

% Individual station images
stn_names = {'CS02', 'CS03', 'CS04', 'CS05', 'CS06', 'CS07', 'CS11', 'CS13', 'CS01', 'CS17', 'CS21', 'CS32'};

if (stncal == 1)
	fprintf (2, '### Run only with calim.parm.minlambda = 5');
	stnmap      = zeros (12,size (acm_t,1), gridparm.uvpad, gridparm.uvpad);
	stncalmap   = zeros (12,size (acm_t,1), gridparm.uvpad, gridparm.uvpad);
	stnuncalmap = zeros (12,size (acm_t,1), gridparm.uvpad, gridparm.uvpad);
	for ind = [1:12]
		stnind   = [(ind-1)*48+1:(ind)*48];
		stnflind = setdiff (stnind, flagant);
		flagant_stn = setdiff ([1:576], stnflind);
		[uloc_x, vloc_x, wloc_x] = gen_flagged_uvloc(uloc, vloc, wloc, flagant_stn);
		[uloc_y, vloc_y, wloc_y] = gen_flagged_uvloc(uloc, vloc, wloc,flagant_stn);
	    for  ts = 1:1 % size (acm_t,1)
	        acm = conj(squeeze(acm_t(ts, stnflind, stnflind,1)));
	        stn(ind,ts)  = pelican_sunAteamsub (conj(squeeze(acm_t(ts,:,:,1))), tmjdsec(ts), fobs, eye(576), flagant_stn,  2, 1,[], [], 'poslocal_afaac12_outer.mat', [], []);
	    
	        [stnmap(ind,ts,:,:), ~, l, m] = fft_imager_sjw_radec (stn(ind,ts).calvis(:), uloc_x(:), vloc_x(:), gridparm, [], [], tmjdsec(ts), fobs, 0);
	        [stnuncalmap(ind,ts,:,:), ~, l, m] = fft_imager_sjw_radec (acm(:), uloc_x(:), vloc_x(:), gridparm, [], [], tmjdsec(ts), fobs, 0);
	    
	        calmat = (stn(ind,ts).gainsol(stnflind)' * stn(ind,ts).gainsol(stnflind)); 
	        [stncalmap(ind,ts,:,:), ~, l, m] = fft_imager_sjw_radec (calmat.*acm, uloc_x(:), vloc_x(:), gridparm, [], [], tmjdsec(ts), fobs, 0);
	    
	    end;
	end;
else
	% Now calibrate the a12 correlations with the station level gains.
	stn_gain_vec = zeros (1, 576);
	for ind = 1:12
	  stn_gain_vec = stn_gain_vec + stn(ind, 1).gainsol;
	end;
	stn_gain_mat = stn_gain_vec' * stn_gain_vec;
	acm_stncal = acm_uncal .* stn_gain_mat;
	
	% Now carry out another level of calibration with the precalibrated vis.
	% For AARTFAAC-6
	fprintf (2, '### Run only with calim.parm.minlambda = 10');
	a6_cal  = pelican_sunAteamsub (acm_uncal, tmjdsec(ts), fobs, eye(576), union([289:576],flagant),  2, 1,[], [], 'poslocal_afaac12_outer.mat', [], []);
	
	a6_stncal  = pelican_sunAteamsub (acm_stncal, tmjdsec(ts), fobs, eye(576), union([289:576],flagant),  2, 1,[], [], 'poslocal_afaac12_outer.mat', [], []);
	
	% For AARTFAAC-12
	a12_cal  = pelican_sunAteamsub (acm_uncal, tmjdsec(ts), fobs, eye(576), flagant,  2, 1,[], [], 'poslocal_afaac12_outer.mat', [], []);
	
	a12_stncal  = pelican_sunAteamsub (acm_stncal, tmjdsec(ts), fobs, eye(576), flagant,  2, 1,[], [], 'poslocal_afaac12_outer.mat', [], []);
	
	
	% Create images
	[a6_cal_img, ~, l, m] = fft_imager_sjw_radec (a6_cal.calvis(:), uloc_a6(:), vloc_a6(:), gridparm, [], [], tmjdsec(1), fobs, 0); 
	
	[a6_stncal_img, ~, l, m] = fft_imager_sjw_radec (a6_stncal.calvis(:), uloc_a6(:), vloc_a6(:), gridparm, [], [], tmjdsec(1), fobs, 0); 
	
	[a12_cal_img, ~, l, m] = fft_imager_sjw_radec (a12_cal.calvis(:), uloc_a12(:), vloc_a12(:), gridparm, [], [], tmjdsec(1), fobs, 0); 
	
	[a12_stncal_img, ~, l, m] = fft_imager_sjw_radec (a12_stncal.calvis(:), uloc_a12(:), vloc_a12(:), gridparm, [], [], tmjdsec(1), fobs, 0); 
end;
% save ('stn_cal_uncal_sol.mat', 'stn', 'stnmap', 'stnuncalmap', 'stncalmap', '-v7.3');

%{
% Test cases:
testnames = {'A6','A6+7', 'A6+8','A6+9','A6+10','A6+11','A6+12','A6+11+12','A6+10+11+12','A6+9+11+12','A6+8+11+12','A6+7+8+9+10','A12'};
flagant_case =  ...
{union([289:576], flagant), ... %          A6
 union([337:576], flagant), ... %          A6 + 7
 union([289:336, 385:576], flagant), ... % A6 + 8
 union([289:384, 433:576], flagant), ... % A6 + 9
 union([289:432, 481:576], flagant), ... % A6 + 10
 union([289:480, 529:576], flagant), ... % A6 + 11
 union([289:528], flagant), ...          % A6 + 12
 union([289:480], flagant), ...          % A6 + 11 + 12
 union([289:432], flagant), ...          % A6 + 10 + 11 + 12
 union([289:384, 433:480], flagant), ... % A6 + 09 + 11 + 12
 union([289:336, 385:480], flagant), ... % A6 + 08 + 11 + 12
 union([481:576], flagant), ... % A6 + 07 + 08 + 09 + 10
 flagant};                     % A12
% A-6 stations only

calmap = zeros (length (flagant_case), gridparm.uvpad, gridparm.uvpad);1
for ind = 1:length (flagant_case)
    tic;
    flagant_stn = flagant_case {ind};
    sol(ind)  = pelican_sunAteamsub (conj(squeeze(acm_t(1,:,:,1))), tmjdsec(1), fobs, eye(576),  flagant_stn, 0, 1,[], [], 'poslocal_afaac12_outer.mat', [], []);
	[uloc_x, vloc_x, wloc_x] = gen_flagged_uvloc(uloc, vloc, wloc, flagant_stn);
	[uloc_y, vloc_y, wloc_y] = gen_flagged_uvloc(uloc, vloc, wloc,flagant_stn);

    used_ants = setdiff ([1:576], flagant_stn);
    acm =  conj(squeeze(acm_t(1,used_ants,used_ants,1)));
    [calmap(ind,:,:), ~, l, m] = fft_imager_sjw_radec (sol(ind).calvis(:), uloc_x(:), vloc_x(:), gridparm, [], [], tmjdsec(1), fobs, 0); 
    toc;
    imagesc (squeeze(calmap(ind,:,:))); title (sprintf ('Case %d', ind));colorbar;
    plot (uloc_x, vloc_x, '.');
end;
%}

