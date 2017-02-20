% Script to read in the components of a standard LOFAR skymodel file, as used 
% in LOFAR calibration.
% Sample header taken from /globaldata/COOKBOOK/Models/A-Team_4.sky:
% # (Name, Type, Patch, Ra, Dec, I, Q, U, V, ReferenceFrequency='7.38000e+07',
% # SpectralIndex='[]', MajorAxis, MinorAxis, Orientation) = format

%, , CygA, 00:00:00, +00.00.00
%CygA_4_1, GAUSSIAN, CygA, 19:59:30.691, +40.43.54.226, 4.263e+03, 0.0, 0.0, 0.0,
%7.38000e+07, [-0.8], 1.82107e+01, 1.16946e+01, 6.
%91353e+01

function [model] = readlofarskymodel (fname, freq)

    fid = fopen (fname);
    dat = [];
    
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
            dat(patch_ind).name = {'CygA'};
        else
            dat(patch_ind).name = dat(patch_ind-1).patch{3}(end);
        end;
        dat(patch_ind).patch = textscan (fid, '%s %s %s %f:%f:%f %s %f %f %f %f %f [%f] %f %f %f','Delimiter', ',', 'CommentStyle', '#');
        patch_ind = patch_ind + 1;
    end;
    fprintf (1, '<-- Found %d patches in skymodel file.\n', patch_ind - 1);
    for ind = 1:patch_ind-1
        fprintf(1, '<--  %s, %d components.\n', dat(ind).name{1}, length(dat(ind).patch{1}) - 1)
    end;

    % Reshape cells to be of the same size
    % NOTE: The minus 1 is because the last patch
    % does not have mismatched sizes, since there is no following patch header.
    for iind = 1:length (dat)-1 
        for jind = 1:12 % Hardcoded 12 due to file structure (although not true
                        % for all skymodel files.
            dat(iind).patch{jind} = dat(iind).patch{jind}(1:end-1);
        end;
    end;

    % Rescale the component fluxes using the spectral index.
    % Add extracted RA/DEC 
    if (~isempty(freq))
	    for ind = 1:patch_ind - 1
            % Flux scaled by spectral index.
            dat(ind).patch{8} = dat(ind).patch{8} ...
                                .*(freq./dat(ind).patch{12}).^dat(ind).patch{13};
            % RA now in radians.
            dat(ind).patch{17} = (dat(ind).patch{4} + dat(ind).patch{5}/60.  ...
                                + dat(ind).patch{6}/3600.) .* (pi/12);

            % Convert DEC to radians
	        for comp = 1:length(dat(ind).patch{1}) % For each individual component
                rem = dat(ind).patch{7}(comp);
                dec = 0;
                % split () and strsplit () do not exist in R2012!
                for ind = 1:3
                    [tok, rem] = strtok (rem, '.');
                    dec = dec + str2double (tok)/60^(ind-1);
                    if ind == 2
                        dec = dec + str2double (rem{1}(2:end))/60^(ind);
                    end
                end;
	            dat(ind).patch{18}(comp) = dec*(pi/180);
            end;
	    end;
    end;
    model = dat;
end
