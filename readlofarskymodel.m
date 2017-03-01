% Script to read in the components of a standard LOFAR skymodel file, as used 
% in LOFAR calibration.
% Sample header taken from /globaldata/COOKBOOK/Models/A-Team_4.sky:
% # (Name, Type, Patch, Ra, Dec, I, Q, U, V, ReferenceFrequency='7.38000e+07',
% # SpectralIndex='[]', MajorAxis, MinorAxis, Orientation) = format

%, , CygA, 00:00:00, +00.00.00
%CygA_4_1, GAUSSIAN, CygA, 19:59:30.691, +40.43.54.226, 4.263e+03, 0.0, 0.0, 0.0,
%7.38000e+07, [-0.8], 1.82107e+01, 1.16946e+01, 6.
%91353e+01

% Arguments:
%  fname    : LOFAR skymodel filename in the above format.
%  freq     : Frequency of operation.
%  res      : Pixel spatial resolution for creating model, arcsec
%  extent   : Pixel extent out to which  the model will be created, arcsec

% Returns:
%  mod   : Cell array containing fields read from the sky model file.
%  rarng : Range of positions corresponding to the cumulative model for each src.
%  decrng: Range of positions corresponding to the cumulative model for each src.
%  modimg: Array of model images generated over the extent and res. specified by
%  the user.
%  

function [mod, rarng, decrng, modimg] = readlofarskymodel (fname, freq, res, extent,  debug)

    fid = fopen (fname);
    mod = [];
    
    if (isempty (extent))
        extent  = 4800; % arcsec, extent of composite image plane patch.
    end;

    if (isempty (res))
        res = 10 ; % arcsec, resolution of composite image plane patch.
    end;

    % Get rid of empty lines and comments.
    while (1)
        lin = fgetl (fid);
        if isempty (lin)
            continue;
        elseif lin(1) == '#'
            continue;
        else
            break;
        end;
    end;

    if (lin(1) == ',') % Have reached the starting of a patch

        patch_ind = 1;
        while (~feof (fid))
            if (patch_ind == 1)
                mod(patch_ind).name = {'CygA'};
            else
                mod(patch_ind).name = mod(patch_ind-1).patch{3}(end);
            end;
            mod(patch_ind).patch = textscan (fid, ...
            '%s %s %s %f:%f:%f %s %f %f %f %f %f [%f] %f %f %f', ...
            'Delimiter', ',', 'CommentStyle', '#');
            patch_ind = patch_ind + 1;
        end;

        fprintf (1, '<-- Found %d patches in skymodel file.\n', patch_ind - 1);
        for ind = 1:patch_ind-1
            fprintf(1, '<--  %s, %d components.\n', mod(ind).name{1},...
                     length(mod(ind).patch{1}) - 1)
        end;
    
        % Reshape cells to be of the same size
        % NOTE: The minus 1 is because the last patch does not have mismatched 
        % sizes, since there is no following patch header.
        for iind = 1:length (mod)-1 
            % Hardcoded 12 due to file structure (although not true
            % for all skymodel files.
            for jind = 1:12 
                mod(iind).patch{jind} = mod(iind).patch{jind}(1:end-1);
            end;
        end;
    
        % Rescale the component fluxes using the spectral index.
        % Add extracted RA/DEC 
        if (~isempty(freq))
    	    for ind = 1:patch_ind - 1
                % Flux scaled by spectral index.
                mod(ind).patch{8} = mod(ind).patch{8} ...
                               .*(freq./mod(ind).patch{12}).^mod(ind).patch{13};
                % RA now in radians.
                mod(ind).patch{17} = (mod(ind).patch{4} ...
                  + mod(ind).patch{5}/60.  + mod(ind).patch{6}/3600.) .* (pi/12);
    
                % Convert DEC to radians
                % For each individual component
    	        for comp = 1:length(mod(ind).patch{1}) 
                    rem = mod(ind).patch{7}(comp);
                    dec = 0;
                    % split () and strsplit () do not exist in R2012!
                    for jind = 1:2
                        [tok, rem] = strtok (rem, '.');
                        dec = dec + str2double (tok)/60^(jind-1);
                        if jind == 2
                            dec = dec + str2double (rem{1}(2:end))/60^(jind);
                        end
                    end;
    	            mod(ind).patch{18}(comp) = dec*(pi/180);
                end;
    	    end;
        end;

    % rarng = {};
    % decrng = {};
    % modimg={};
    % Now create a spatial grid with the specified resolution for holding all
    % components of a patch cumulatively.
    for patch = 1:length(mod)
        fprintf (1, '<-- Working on modeling response to patch %s (%.0f components).\n', ...
                mod(patch).name{1}, length(mod(patch).patch{1}));
        meanra  = mean (mod(patch).patch{17}); % In rad.
        meandec = mean (mod(patch).patch{18}); % In rad.
        extent_rad = (extent/2)*(pi/180/3600); % Convert extent to rad.
        extent_pix = extent/res;
        res_rad = ((res/3600)*pi/180);
        rarng {patch} = meanra  - extent_rad:res_rad:meanra  + extent_rad;
        decrng {patch}= meandec - extent_rad:res_rad:meandec + extent_rad;
        modimg {patch} = zeros (length(rarng{patch}));

        % nrows = size (pos, 1);
        % ncols = size (pos, 2);
        % [l,m] = radectolm (meanra, meandec, tobs_jd, rodata.lon*180/pi, ...
        %                     rodata.lat*180/pi, false);
        fprintf (1, '<-- Patch mean ra/dec:%.4f, %.4f\n',meanra, meandec);

        tot_flux = 0;
        for jind = 1:length(mod(patch).patch{1}) % Get number of components.
            comp = mod(patch).patch;
            tmp = comp{2}(jind); % Can't use comp(2){jind} directly! :/

            % Offset of this component from the patch mean, in pixel units.
            offra    = int32(((comp{17}(jind) - meanra)*180*3600/pi)/res);
            offdec   = int32(((comp{18}(jind) - meandec)*180*3600/pi)/res);

            % Generate the positions of all the points making up the source.
            if (strcmp (tmp{1},'POINT'))
                modimg{patch} (offdec + int32(extent_pix/2), ...
                       offra  + int32(extent_pix/2)) = ...
                modimg{patch} (offdec + int32(extent_pix/2), ...
                       offra  + int32(extent_pix/2)) + comp{8}(jind);

            elseif (strcmp (tmp{1},'GAUSSIAN'))
                rasig  = (comp{14}(jind)/3600)*(pi/180);
                decsig = (comp{15}(jind)/3600)*(pi/180);
                if (rasig == 0 && decsig == 0)
                    rasig  = (1/3600)*(pi/180); % Set to 1 arcsec (arbit)
                    decsig = (1/3600)*(pi/180); % Set to 1 arcsec (arbit)
                elseif (rasig == 0)
                    rasig = decsig;
                elseif (decsig == 0)
                    decsig = rasig;
                end;
                
                [rapos, decpos, flux] = gengaussiansrc (comp{17}(jind), ...
                   comp{18}(jind), true, comp{8}(jind), rasig, decsig, ...
                   comp{16}(jind), (res/3600)*(pi/180));
                % pos = radectoITRF ((rapos(:)), (decpos(:)), true, tobs_jd);

                decstart = offdec-int32(size (flux, 1)/2) + int32(extent_pix/2); 
                rastart  = offra -int32(size (flux, 2)/2) + int32(extent_pix/2);
    
                modimg{patch} (decstart:decstart + size (flux,1)-1, ...
                       rastart :rastart  + size (flux,2)-1) = ...
                      modimg{patch} (decstart:decstart + size (flux,1)-1, ...
                             rastart :rastart  + size (flux,2)-1) + flux;
            end;


            tot_flux = tot_flux + comp{8}(jind);
            
            % simsky_up = pos * rodata.normal > 0;
            % if (sum (simsky_up) == 0)
            %     fprintf (2, 'Not visible.');
            %    continue;
            % end;

	        % simsky_A = exp(-(2*pi*1i*freq/rodata.C)*(rodata.posITRF_fl*pos.'));
	        % simsky_acc = simsky_acc +  diag(antgains) * simsky_A*diag(flux(simsky_up))* ...
		    % 		simsky_A' * diag (antgains)' + sigman;
	        % simsky_acc = simsky_acc +  simsky_A*diag(flux(simsky_up))* ...
		% 		simsky_A';
            
            % fprintf (1, '%d grid(%d,%d)   ', jind, size(flux,1), size(flux,2));
            fprintf (1, '%d ', jind);
        end;

        if (debug > 4)
            figure ;
            imagesc (rarng{patch}, decrng{patch}, modimg{patch}); colorbar;
            title (sprintf ('Model of %s', mod(patch).name{1}));
        end;
        fprintf (1, '..Done.\n Total flux: %f.\n', tot_flux);
    end;
end;
end

% Convert a gaussian model source to positions and fluxes 
% Arguments:
%  ra_cen : RA of gaussian peak, rad.
% dec_cen : Dec of gaussian peak, rad.
% epoch   : True => J2000
% tobs_jd : Time of observation in JD
% pk_flux : Peak flux of component, Jy.
% sigr/d  : Gaussian width in RA and DEC, in radians.
% pa      : Position angle in deg, in radians.
% res     : Resolution at which the Gaussian is to be constructed, in radians.
%
% Returns:
%   ragrid,decgrid: The positions of all parts of the component in radians    
function [ragrid, decgrid, flux] = gengaussiansrc (ra_cen, dec_cen, epoch, pk_flux, sigr, sigd, pa, res)

    extent = 3; % Extent of gaussian in sigma units.
    % Create a grid of points with resolution res radians, to cover the gaussian
    rarng = ra_cen + [-extent*sigr:res:extent*sigr];
    decrng = dec_cen + [-extent*sigd:res:extent*sigd];
    [ragrid, decgrid] = meshgrid (rarng, decrng);
    flux = pk_flux * exp (-( (ragrid-ra_cen).^2/(2*sigr^2) + (decgrid-dec_cen).^2/(2*sigd^2) ));
    % return ragrid, decgrid, flux;
end
