% Script to generate an image after integrating in the visibility domain
cd '/dop312_0/prasad/GPU_CORR_DAT/afaac-6/8sb_29Jan16/'
fnames = {'SB295.vis', 'SB296.vis', 'SB297.vis', 'SB298.vis', 'SB299.vis', 'SB300.vis', 'SB301.vis', 'SB302.vis'};
obs.sub = [295:302];
nsub = 8;
obs.flagant_x = [18, 84, 142, 168, 262];
obs.flagant_y = [18, 84, 142, 168, 262];
obs.cal = 1;
obs.chans = 1:63;
obs.chanwidth = 195312.5/64;
obs.nelem = 288
obs.uvflag_x = eye(obs.nelem);
obs.uvflag_y = eye(obs.nelem);
obs.deb = 0;
obs.ptSun = 1;
obs.posfilename = 'poslocal_outer.mat';
fout = './';

for i = [1:nsub]
    obs.freq = obs.sub(i)*195312.5;
    fprintf (1, '<-- Creating VisRec object for file %s with freq %f.\n', fnames{i}, obs.freq);
    sbrecobj (i) = VisRec(fnames{i}, obs); 
end;


% Calibrate all subbands
    acm_x = zeros (nsub, obs.nelem, obs.nelem);
    acm_y = zeros (nsub, obs.nelem, obs.nelem);
    t1 = triu (ones (obs.nelem));

    for sb = [1:nsub]
        sbrecobj(sb).readRec ([1,0,0,1], obs.chans);
        acm_tmp = zeros (obs.nelem);
        acm_tmp (t1 == 1) = sbrecobj(sb).xx;
        acm_tmp = acm_tmp + acm_tmp';
        acm_tmp (eye(obs.nelem) == 1) = real(diag(acm_tmp));
        acm_x(sb,:,:) = conj(acm_tmp);

        fprintf (2, '\n<-- Calibrating XX for subband %d.\n', sb);
        sol_x(sb) = pelican_sunAteamsub (squeeze(acm_x(sb,:,:)), sbrecobj(sb).trecstart, sbrecobj(sb).freq, obs.uvflag_x, obs.flagant_x, obs.deb, obs.ptSun, [], [], obs.posfilename, [], []);

        acm_tmp = zeros (obs.nelem);
        acm_tmp (t1 == 1) = sbrecobj(sb).yy;
        acm_tmp = acm_tmp + acm_tmp';
        acm_tmp (eye(obs.nelem) == 1) = real(diag(acm_tmp));
        acm_y(sb,:,:) = conj(acm_tmp);
        
        fprintf (2, '\n<-- Calibrating YY for subband %d.\n', sb);
        sol_y(sb) = pelican_sunAteamsub (squeeze(acm_y(sb,:,:)), sbrecobj(sb).trecstart, sbrecobj(sb).freq, obs.uvflag_y, obs.flagant_y, obs.deb, obs.ptSun, [], [], obs.posfilename, [], []);
    end;

    % Local horizon based coordinates of array in ITRF
    load (obs.posfilename, 'posITRF', 'poslocal'); 

    % Generate uv coordinates in local horizon coord. system, needed for imaging
    uloc = meshgrid (poslocal(:,1)) - meshgrid (poslocal (:,1)).';
    vloc = meshgrid (poslocal(:,2)) - meshgrid (poslocal (:,2)).';

    [uloc_x, vloc_x] = gen_flagged_uvloc (uloc, vloc, obs.flagant_x);
    [uloc_y, vloc_y] = gen_flagged_uvloc (uloc, vloc, obs.flagant_y);
    obs.gridparm.type = 'pillbox';
    obs.gridparm.lim  = 0;

    % Control imaged field of view, independently of frequency.
    % obs.gridparm.duv  = (0.5/180)*obs.fov;
    obs.gridparm.duv  = 0.5;

    obs.gridparm.Nuv  = 1024;    % size of gridded visibility matrix
    obs.gridparm.uvpad= 1024;    % specifies if any padding needs to be added
    obs.gridparm.fft  = 1;

% Add up all the subbands in visibility plane and image
for sb = 1:nsub
    integvis_x = zeros (size (sol_x(1).calvis));
    integvis_y = zeros (size (sol_y(1).calvis));
    for t = 1:sb
        integvis_x = integvis_x + sol_x(t).calvis;
        integvis_y = integvis_y + sol_y(t).calvis;
    end;
    integvis_x = integvis_x ./ sb;
    integvis_y = integvis_x ./ sb;
    [~, integmap_x, ~, l, m] = ... 
        fft_imager_sjw_radec (integvis_x(:), uloc_x(:), vloc_x(:), ... 
            obs.gridparm, [], [], sbrecobj(sb).trecstart, mean(obs.sub(1:sb))*195312.5, 0);
    [~, integmap_y, ~, l, m] = ... 
        fft_imager_sjw_radec (integvis_y(:), uloc_y(:), vloc_y(:), ... 
            obs.gridparm, [], [], sbrecobj(sb).trecstart, mean(obs.sub(1:sb))*195312.5, 0);


    % Generate FITS filename for this final image
    img.tobs       = sbrecobj(1).trecstart;
    img.pix2laxis = size (integmap, 1);
    img.pix2maxis = size (integmap, 2);
    img.freq       = mean (obs.sub*195312.5);
    img.df        = sb*length (obs.chans)*obs.chanwidth;
    img.map       = (integmap_x + integmap_y)/2;;
    
    integmapname = sprintf ('%s/Sb%3d-%3d_R%02d-%02d_T%s.fits', fout, obs.sub(1), obs.sub(sb), obs.chans(1), obs.chans(end), datestr (mjdsec2datenum (sbrecobj(1).trecstart), 'dd-mm-yyyy_HH-MM-SS')); 
    wrimg2fits (img, integmapname);
end;
