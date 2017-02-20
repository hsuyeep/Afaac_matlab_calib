% Script to generate the power patterns of individual dipoles from the
% simulation of the 6 station system by Arts
% pep/28Jul15

% NOTE: Moved computation to SimElementBeam.m. Keeping for historical significance.
% pep/17Nov15. 

% Load the simulation outputs. See genafaacbeam.m for a description of the
% data products.
load ('LBA_beam/arts_sim/LBA_core.mat');      % 5MHz grid, 10-90MHz
load ('LBA_beam/arts_sim/LBA_core_fine.mat'); % 1MHz grid,Load the simulation products
% load ('LBA_beam/arts_sim/HBA_station.mat');   % HBA, currently untested.

fprintf (1, '<-- Loaded simulation products...\n');
% ant = [1, 1151];
% array = 'lba_outer';
% 
% % Generate a selection of dipoles based on the station array configuration
% switch array
% 		case 'lba_inner'    % The first 48 of 96 dipoles per station
% 			sel_x = [1:2:576]; % The X-dipole of the LBA_INNER antennas for the 6 stations
%             sel_y = sel_x + 1;
% 
% 		case 'lba_outer'
% 			sel_x = [577:2:1152]; % The X-dipole of the LBA_OUTER antennas for the 6 stations
%             sel_y = sel_x + 1;			
% end;
% ant = sel_x;    
 
% Generate output voltages from LNAs, per dipole, per theta/phi.
	% Nelem = length (sel); % Number of selected antenna elements.
	Nelem = 1152;
	Glna  = 1; % Define Voltage gains of the LNAs, arbit. as unity.
	Vout  = zeros(Nelem, 2, length(Theta), length(Phi), length(Freq));
	
    for NF=1:length(Freq)
        % Zlna taken from the Master's thesis of Maria Krause, pg. 27
        % (Schematic) and pg. 47 for estimated values of R and C.
		Zlna = (700./(1+i*2*pi*Freq(NF)*700*15e-12))*eye (Nelem); % Define impedances of the LNAs.		
		B=Glna*Zlna/(Zant(:,:,NF)+Zlna); % Is it correct to ignore the other cross terms?
		for NT=1:length(Theta)
			for NP=1:length(Phi)
	       		Vout(:,:,NT,NP,NF)=B*g_z(:,:,NT,NP,NF);				
			end
        end		
    end
    fprintf (1, '<-- Generated output voltages per dipole, per E-field component.\n');
    % Generate Gain matrix (see header comment in genafaacbeam.m) per
    % antenna, for all antennas and both pols.
    G = zeros (length(Theta), length(Phi), length(Freq), Nelem/2, 2, 2);
    
    for NF=1:length(Freq) 
        for NT=1:length(Theta)
            for NP=1:length(Phi)
                for ant=1:2:Nelem % Because we use up both pols of the ant.
                G(NT, NP, NF, (ant+1)/2, :, :) = [squeeze(Vout(ant  ,1,NT,NP,NF)) squeeze(Vout(ant  ,2,NT,NP,NF));
                  squeeze(Vout(ant+1,1,NT,NP,NF)) squeeze(Vout(ant+1,2,NT,NP,NF))];
                end;
            end;
        end;
    end;
    fprintf (1, '<-- Calculated full pol. Jones gain matrix per antenna.\n');
    
    % Generate Coherency matrix and stokes-I beam for all antennas
    CM = zeros (size (G));
    stI = zeros (Nelem/2, length (Theta), length(Phi), length(Freq));    
    for NF = 1:length(Freq)
        for NT=1:length(Theta)
            for NP=1:length(Phi)
                for ant = 1:Nelem/2
                    tmp = squeeze(G(NT,NP,NF,ant,:,:)); % 2x2 Jones matrix
                    CM(NT,NP,NF,ant,:,:) = tmp*0.5*eye(2)*tmp';
                    % below is already a power due to adding up xx and yy
                    % of the CM.
                    stI (ant,NT, NP, NF) = CM(NT,NP,NF,ant,1,1) + CM(NT,NP,NF,ant,2,2); 
                end;
            end;
        end;
    end;
    fprintf (1, '<-- Generated coherency matrix and stokes-I power beam.\n');
    
    % Generate the interpolated beam on the l,m plane
    l = [-1:0.01:1]; m = l;
    
    [lgrid, mgrid] = meshgrid(l, m);
	dist = sqrt(lgrid.^2 + mgrid.^2);
	thetai = asin(dist);
	thetai(dist >= 1) = 0;
	thetai = thetai * 180/pi; % Convert to degrees.
	phii = mod(atan2(lgrid, mgrid), 2 * pi)*180/pi; % Convert to degrees.

	% Interpolate the simulated beam to the frequency and l,m positions needed 
	% by the user, for all 576 antenna
	beam = zeros(Nelem/2, length(l), length(m), length(Freq));
	for idx = 1:length(Freq)
        for dip = 1:Nelem/2
            beam(dip,:, :, idx) = interp3 (Phi, Theta, Freq, squeeze(stI(dip,:,:,:)), phii, thetai, Freq(idx) * ones(size(phii))) .* (dist < 1);
        end;
    end
    fprintf (1, '<-- Generated interpolated beam.\n');
    
%     % Plot the beam for the first antenna
%      for ind = 1:length(Freq)
%         imagesc (squeeze(abs(beam(1,:,:,ind)))./squeeze(max(max(abs(beam(1,:,:,ind)))))); colorbar;
%         title (sprintf ('Stokes-I beam for Antenna %d, freq %dMHz', 1, Freq(ind)/1e6))
%         pause(1);
%     end;
%     
    % Plot the average beam for all stations, LBA_OUTER.
    % NOTE: We can separate out the beams because the mutual coupling
    % exists as long as the antennas are in the vicinity, and does not
    % depend on them being turned on or not.
    figure ('PaperType', 'a4');
    set(0, 'DefaultTextFontSize', 12);
    set(0, 'DefaultAxesFontSize', 12);
    
    tmp1 = squeeze(mean (reshape (beam (  1: 288,:,:,:), [48, 6, 201, 201, 17]), 1));
    tmp2 = squeeze(mean (reshape (beam (289: 576,:,:,:), [48, 6, 201, 201, 17]), 1));
    beam_inner = zeros (size (tmp1));
    beam_outer = beam_inner;
    for ind=1:6
        for iind = 1:length(Freq)
            beam_inner (ind, :, :, iind) = tmp1(ind,:,:,iind)./max(max(squeeze(abs(tmp1(ind,:,:,iind)))));
            beam_outer (ind, :, :, iind) = tmp2(ind,:,:,iind)./max(max(squeeze(abs(tmp2(ind,:,:,iind)))));
        end;
    end;
    clear tmp1, tmp2;
    
        
    station_names = {'CS02', 'CS03', 'CS04', 'CS05', 'CS06', 'CS07'};
    for idx = 1:length(Freq)        
        for station = 1:6                 
            subplot (6,2,2*station-1);
            b1 = abs(squeeze(beam_outer(station, :,:,idx)));
            imagesc (l,l,b1); colorbar;
            
            c = contourc (l, l, b1, [0.75, 0.5, 0.25]);
            hold on;
            c1 = c(2,1); c2 = c(2,c1+2);
            plot (c(1,2:c1+1), c(2,2:c1+1), '-w');
            plot (c(1, c1+3:c1+c2+2), c(2,c1+3:c1+c2+2), '-k');
            plot (c(1, c1+c2+4:end), c(2,c1+c2+4:end), '-c');
            %axis equal
            % axis tight
            set (gca, 'YDir', 'Normal'); % To match orientation with station images
            set (gca, 'XDir', 'Reverse'); % To match orientation with station images
            text (-0.70, -0.75, sprintf ('%s', station_names{station}), 'Color', [1,1,1]);            
            hold off;
            
            
            
            subplot (6,2,2*station);
            b2 = abs(squeeze(beam_inner(station, :,:,idx)));
            imagesc (l,l,b2); colorbar;
            c = contourc (l, l, b2, [0.75 0.5, 0.25]);
            hold on;
            c1 = c(2,1); c2 = c(2,c1+2);
            plot (c(1,2:c1+1), c(2,2:c1+1), '-w');
            plot (c(1, c1+3:c1+c2+2), c(2,c1+3:c1+c2+2), '-k');
            plot (c(1, c1+c2+4:end), c(2,c1+c2+4:end), '-c');
            % axis equal
            % axis tight
            set (gca, 'YDir', 'Normal'); % To match orientation with station images
            set (gca, 'XDir', 'Reverse'); % To match orientation with station images
            text (-0.70, -0.75, sprintf ('%s', station_names{station}), 'Color', [1,1,1]);            
            hold off;
        end;               
        mtit (sprintf ('Stokes-I power patterns, %0.2f MHz, LBA-OUTER (left), LBA-INNER(right)', Freq(idx)/1e6), 'xoff', -0.1, 'yoff', 0.025);
        print ('-depsc', '-r200', sprintf ('sterp_stations_stokesI_%0.2fMHz.png', Freq(idx)/1e6));
        % saveas (gcf, sprintf ('sterp_stations_stokesI_%0.2fMHz.png', Freq(idx)/1e6));
        pause;     
        clf;
    end;        
    fprintf (1, '<-- Plotted mean primary beam per station.\n');
    
    % Generate the spectral response per station at the zenith and various elevations.    
    % The l/m coordinates of the .75, 0.5, 0.25 contour at 60
    % MHz/LBA_outer, converted to pixel coordinates.
    sample_beam = [0.01, 0.5, 0.7, 0.9];
    sample_pix = int32 (length(l)/2) + int32(sample_beam*length(l)/2);
    
    beam_inner_freq = zeros (6, length (Freq), length(sample_pix));
    beam_outer_freq = zeros (6, length (Freq), length(sample_pix));
    
    for station = 1:6
        for ind = 1:length(Freq)
            beam_inner_freq (station, ind, :) = diag(squeeze(abs(beam_inner(station, sample_pix, sample_pix, ind))));
            beam_outer_freq (station, ind, :) = diag(squeeze(abs(beam_outer(station, sample_pix, sample_pix, ind))));
        end;
    end;
    
    clf;
    for ind = 1:6
        subplot (6, 2, 2*ind-1);
        plot (Freq/1e6, squeeze(beam_inner_freq(ind, :, :)), '.-')
        xlabel ('Freq(MHz)'); ylabel ('Normalized linear gain');
        title (sprintf('%s:LBA-INNER', station_names{ind}));
        
        subplot (6, 2, 2*ind)
        plot (Freq/1e6, squeeze(beam_outer_freq(ind, :, :)), '.-')
        xlabel ('Freq(MHz)'); ylabel ('Normalized linear gain');
        title (sprintf('%s:LBA-OUTER', station_names{ind}));        
    end;
    mtit (sprintf ('Beam pattern spectral response: %s FoV', num2str(sample_beam))); 
    saveas (gcf, 'beam_spectral_arrayconfig_response.png');
    
    
    % Generate a movie of the stokes-I beam for each antenna, LBA_INNER
    fig1 = figure;
    freqind = 4;
    movname = strcat (sprintf ('LBA_INNER_stokesI_beam_pattern_%fMHz', Freq(freqind)), '.avi');
    vidObj = VideoWriter (movname, 'Motion JPEG AVI');
    vidObj.FrameRate = 5;
    open (vidObj);

    dipsel = [1:2:576];
    for station = 1:6
        for  ind = 1:48 
            clf
            ant = (station-1)*48 + ind;
            stat = (station-1)*48+1;
            subplot (1,3,1);
            plot (Positions(1, dipsel), Positions(2, dipsel), '.');
            hold on;
            plot (Positions(1, 2*ant), Positions(2, 2*ant), 'r*');
            hold off;
            subplot (1,3,2); tmp = squeeze(abs(beam(ant,:,:,freqind))); imagesc (l,l,tmp./max(max(tmp))); colorbar;
            title (sprintf ('CS00%d INNER:Ant %d, %.2f MHz', station+1, ind, Freq(freqind)/1e6));
            
            subplot (1,3,3); tmp = squeeze(mean(squeeze(abs(beam(stat:ant,:,:,freqind))), 1)); imagesc (l,l,tmp./max(max(tmp))); colorbar;
            title (sprintf ('CS00%d INNER Mean beam, %.2f MHz', station+1, Freq(freqind)/1e6));
            pause(0.1);
            currFrame = getframe (fig1);
            writeVideo (vidObj, currFrame);
        end;
    end;
    close(vidObj);
    
    % LBA_OUTER now.
    clf;
    freqind = 4;
    movname = strcat (sprintf ('LBA_OUTER_stokesI_beam_pattern_%fMHz', Freq(freqind)), '.avi');
    vidObj = VideoWriter (movname, 'Motion JPEG AVI');
    vidObj.FrameRate = 5;
    open (vidObj);

    dipsel = [577:2:1152];
  
    for station = 1:6
        for  ind = 1:48 
            clf
            ant = 288 + (station-1)*48 + ind;
            stat = 288 + (station-1)*48+1;
            subplot (1,3,1);
            plot (Positions(1, dipsel), Positions(2, dipsel), '.');
            hold on;
            plot (Positions(1, 2*ant), Positions(2, 2*ant), 'r*');
            hold off;
            subplot (1,3,2); tmp = squeeze(abs(beam(ant,:,:,freqind))); imagesc (l,l,tmp./max(max(tmp))); colorbar;
            title (sprintf ('CS00%d OUTER:Ant %d, %.2f MHz', station+1, ind, Freq(freqind)/1e6));
            
            subplot (1,3,3); tmp = squeeze(mean(squeeze(abs(beam(stat:ant,:,:,freqind))), 1)); imagesc (l,l,tmp./max(max(tmp))); colorbar;
            title (sprintf ('CS00%d OUTER Mean beam, %.2f MHz', station+1, Freq(freqind)/1e6));
            pause(0.1);
            currFrame = getframe (fig1);
            writeVideo (vidObj, currFrame);
        end;
    end;      
    close(vidObj);
    
    % Generate the average beam pattern over all 6 stations. Exclude CS004
    % for the 20 and 40 MHz power patterns (or use median?)    
    aartfaac6_inner_beam = squeeze(mean (beam_inner, 1));
    aartfaac6_outer_beam = squeeze(mean (beam_outer, 1));
 
    aartfaac6_inner_beam (:,:,4) = squeeze(mean (beam_inner ([1,2,4,5,6],:,:, 4)));
    aartfaac6_inner_beam (:,:,7) = squeeze(mean (beam_inner ([1,2,4,5,6],:,:, 7)));
 
    aartfaac6_outer_beam (:,:,4) = squeeze(mean (beam_outer ([1,2,4,5,6],:,:, 4)));
    aartfaac6_outer_beam (:,:,7) = squeeze(mean (beam_outer ([1,2,4,5,6],:,:, 7)));
    % Ignore CS004 for the freq indices 4(=25MHz) and 7(=40MHz)

   
    for iidx = 1:2
        clf
       for idx = 1:9
            ind = (iidx-1)*9 + idx;
            if (ind > length(Freq))
                break;
            end;
            subplot (3,3,idx);
            b1 = abs(aartfaac6_outer_beam(:,:,ind));
            imagesc (l,l, b1); colorbar;
            
            c = contourc (l, l, b1, [0.75, 0.5, 0.25]);
            hold on;
            c1 = c(2,1); c2 = c(2,c1+2);
            plot (c(1,2:c1+1), c(2,2:c1+1), '-w');
            plot (c(1, c1+3:c1+c2+2), c(2,c1+3:c1+c2+2), '-k');
            plot (c(1, c1+c2+4:end), c(2,c1+c2+4:end), '-c');            
            set (gca, 'YDir', 'Normal'); % To match orientation with station images
            set (gca, 'XDir', 'Reverse'); % To match orientation with station images             
            title (sprintf ('%.2fMHz', Freq(ind)/1e6));
            hold off;                                               
        end;
        mtit (sprintf('AARTFAAC-6 Stokes-I power patterns, LBA-OUTER'), 'xoff', -0.1, 'yoff', 0.025);
        saveas (gcf, sprintf('aartfaac6_outer_beam_%d.png', iidx));     
    end;
        
    % plot the spectral variation of the power contours'
    col = {'b', 'g',' r','c', 'm', 'y', 'k', 'w'};
    leginfo = {};
    for pind = [0.75, 0.5, 0.25]
        
        figure;        
        for ind = 1: length(Freq)/2
            b1 = abs(aartfaac6_outer_beam (:,:,2*ind));
            c = contourc(l, l, b1, [pind, pind]);
            fprintf (1, '%d ', c(2,1));
%             subplot (1,2,1);
%             plot (c(1,2:end));
%             subplot (1,2,2);
%             plot (c(2,2:end));
              plot (c(1,2:end), c(2,2:end), char (col{ind}));
            hold on;            
            leginfo{ind} = [num2str(Freq(ind)/1e6) 'MHz']; %sprintf ('%.1MHz', Freq(ind)/1e6) ];           
        end;
        legend (leginfo);
        title (sprintf ('LBA-OUTER, Gain level %f.\n', pind));
        grid on;
        % legend (sprintf ('%.1fMHz,', Freq(1:2:end)/1e6));
        fprintf (1, '\n');
    end;
