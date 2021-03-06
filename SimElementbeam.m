classdef SimElementbeam < handle
    % Class to model the element beam from a simulation.
    % For HBA, this entails beamforming, while for LBA, no beamforming is
    % carried out. For both situations, an open circuit voltage is
    % available (directly from LBA, and after beamforming from the HBA).
    % 
    % - Allows beamforming in a chosen direction.
    % - Generates the tile beam response for the chosen direction
    % - Allows generation of the tile response while tracking the chosen sky direction. 
    % pep/06Aug15
    % Usage:
    %  obj = SimElementbeam (1, 'lbainner', './', []);
    %  obj.createAvgBeam(1, [-1:0.01:1], [-1:0.01:1], []);                         
    %  obj.genPlots ();
    properties
        lba         = 0;  % Bool to indicate whether an LBA or HBA simulation
        simul_file  = []; % File name of file containing the sim. products.
        simloaded   = 0;  % Flag recording whether the simulation has been already loaded.
        arrayconfig = []; % Stores the chosen array configuration.
        dipsel      = []; % Dipole selection subset from the array, based on the arrayconfig.
        nelem       = []; % Number of dipole elements in the array.
        nfreq       = []; % Number of freq elements in the normal dipole.
        
        % Simulation specific datastructure
        glna        = 0;  % The gain of the LNA element.
        zlna        = 0;  % The impedance seen by the antenna LNA;
        vout        = []; % Placeholder for variables with LNA output voltage.
        gain        = [];
        g_z  = [];
        Zant = [];
        Freq = [];
        Phi  = [];
        Theta= [];
        stokes = [];
        usrfreq = [];
        Positions = [];
        radecstokesbeam = [];
        lmstokesbeam = [];
        station_names = {'CS02', 'CS03', 'CS04', 'CS05', 'CS06', 'CS07'};
        fsaveprefix = [];
        lmstnavgbeam = [];
        lmaartfaacbeam = [];
        pointing = []; % Pointing direction for HBA analog beamforming, azi/el, rad        
        l = [];
        m = [];
        fhdl = [];
    end
    
    methods
        
        % Constructor
        % Arguments:
        %  lba : Bool, 1 for an LBA observation, 0 for an HBA observation
        % arrayconfig: one of LBA_OUTER, LBA_INNER, HBA_0/1
        %  fsaveprefix: The prefix with which to save generated plots.
        % pointing : Pointing vector direction for HBA analog beam forming
        % simulation (azi/el, rad)
        function object = SimElementbeam (lba, arrayconfig, fsaveprefix, pointing)
            addpath ('~/WORK/Matlab');
            object.lba = lba;
            object.arrayconfig = arrayconfig;
            assert (strcmpi(arrayconfig, 'LBAOUTER') || strcmpi(arrayconfig,'LBAINNER') || strcmpi(arrayconfig,'HBA0') || strcmpi(arrayconfig, 'HBA1'));
            
            if (object.lba == 1)
                object.simul_file = '/home/prasad/WORK/AARTFAAC/Afaac_matlab_calib/LBA_beam/arts_sim/LBA_core.mat';                
            else
                object.simul_file = '/home/prasad/WORK/AARTFAAC/Afaac_matlab_calib/LBA_beam/arts_sim/HBA_station.mat';                
            end;
            fprintf (2, '<-- Loading simulation products for %s simulations from %s.\n', object.arrayconfig, object.simul_file);
            load (object.simul_file);
            object.g_z  = g_z;
            object.Zant = Zant;
            object.Freq = Freq;
            object.Phi  = Phi;
            object.Theta= Theta;
            object.Positions = Positions; 
            object.nelem = length (object.g_z);
            object.simloaded = 1;
            object.fsaveprefix = fsaveprefix;
            if (object.fsaveprefix(end) ~= '/')
                object.fsaveprefix = [object.fsaveprefix '/'];
            end;
            
            
            if  (object.lba == 1 && strcmpi(object.arrayconfig,'LBAINNER'))
                object.dipsel = 1:576;            % 288 dual pol antennas
            elseif (object.lba == 1 && strcmpi(object.arrayconfig,'LBAOUTER'))
                object.dipsel = 577:1152;         % 288 dual pol antennas
            elseif (object.lba == 0 && strcmpi(object.arrayconfig, 'HBA'))
                object.dipsel= 1:768;             % 16 dual pol ants/tilex24tiles
            else
                fprintf (2, '### Unknown array configuration! Quitting...\n');
                return;
            end;
        
            object.nfreq = length (object.Freq);
            % Create datastructures
            % ALERT! The datastructures, if initialized in the member
            % functions, do not hold their data past the member function
            % call, even if they are declared here, e.g. object.vout = []!
            object.vout    = zeros(int32(object.nelem/2), 2, length(object.Theta), length(object.Phi), object.nfreq);
            object.gain    = zeros (length(object.Theta), length(object.Phi), object.nfreq, object.nelem/4, 2, 2);           
            object.radecstokesbeam = zeros (object.nelem/4, length (object.Theta), length(object.Phi), object.nfreq, 4); 
            object.lmstokesbeam = []; %zeros(object.nelem/4, length(l), length(l), length(freq), 4);
            object.glna = 1;
            object.zlna = 1;
            if (isempty (pointing))
                pointing = [0, pi/2];
            end;
        end;
        
       
%%%%%%%%%%%% Beam generation related helper functions %%%%%%%%%%%%        
        % Generate the expected output voltage of the LNA, given the input
        % open circuit voltage and LNA impedance at the simulated
        % frequencies, at all stokes.
        % Arguments:
        %  None
        % Returns:
        %  None
        function genOutputVoltage (object)
            assert (object.simloaded == 1);
            if (isempty(object.vout) == 1)
                object.vout    = zeros(int32(object.nelem/2), 2, length(object.Theta), length(object.Phi), object.nfreq);
            end;
            for NF=1:object.nfreq
                if (strcmpi(object.arrayconfig,'LBAOUTER') || strcmpi(object.arrayconfig,'LBAINNER'));
                    % Zlna taken from the Master's thesis of Maria Krause, pg. 27
                    % (Schematic) and pg. 47 for estimated values of R and
                    % C.
                    object.zlna = (700./(1+1i*2*pi*object.Freq(NF)*700*15e-12))*eye (object.nelem); % Define impedances of the LNAs.		                                    
                    % object.zlna =(36+1./(1i.*2.*pi.*object.Freq(NF).*4.1e-12)).*eye(object.nelem);% 36 + (1./(1i*2*pi*object.Freq(NF)*4.1e-12)) * eye (object.nelem); % R & C in series, values from LBA amplifier schematics.
                    % object.zlna = 0 + (1./1i*2*pi*object.Freq(NF)*4.1e-12) * eye (object.nelem); % R & C in series, values from LBA amplifier schematics.
                end;
                B = (object.glna*object.zlna(object.dipsel,object.dipsel)/(object.Zant(object.dipsel,object.dipsel, NF)+object.zlna(object.dipsel,object.dipsel))); % Is it correct to ignore the other cross terms?                
                for NT=1:length(object.Theta)
                    for NP=1:length(object.Phi)
                        object.vout(:,:,NT,NP,NF)=B*object.g_z(object.dipsel,:,NT,NP,NF);				
                    end
                end		
            end
            fprintf (1, '<-- Generated output voltages per dipole, per E-field component.\n');
        end;
        
        % Function to generate Jones matrix based on antenna selection.
        % at the simulated frequencies.
        % Arguments:
        %  None
        % Returns:
        %  None
        function genJonesMat (object) 
            assert (object.simloaded == 1);            
            assert (isempty(object.vout) == 0);
            if (isempty(object.vout))
                % 2x2 Jones matrix which incorporates both pols and E
                % components.
                object.gain    = zeros (length(object.Theta), length(object.Phi), object.nfreq, object.nelem/4, 2, 2);
            end;
            for NF=1:object.nfreq
                for NT=1:length(object.Theta)
                    for NP=1:length(object.Phi)
                        % Because we use up both pols of the ant, and only half the antennas are 
                        % used in any LBA_INNER/OUTER configuration.
                        for ant=1:2:object.nelem/2 
                            object.gain(NT, NP, NF, (ant+1)/2, :, :) = [squeeze(object.vout(ant  ,1,NT,NP,NF)) squeeze(object.vout(ant  ,2,NT,NP,NF));  % [g_xth g_xph;
                                                                        squeeze(object.vout(ant+1,1,NT,NP,NF)) squeeze(object.vout(ant+1,2,NT,NP,NF))]; %  g_yth g_yph];                           
                        end;
                    end;
                end;
            end;
        
            fprintf (1, '<-- Calculated full pol. Jones gain matrix per antenna.\n');            
        end;              
        
        
        % Function to generate actual beams, for LBA antennas, at the
        % simulated frequencies.
        % Arguments:
        %  stokes: Cell array of which stokes beams are required
        % Returns:
        %  None: Generated outputs are internally stored, in the same order
        %  as specified in the stokes array.
        function genRaDecStokesBeam (object, stokes)
            assert (object.simloaded == 1);      
            assert (isempty(object.gain) == 0)
            if (isempty(stokes))
                stokes = 'I'; % Generate stokes-I by default.
            end;
            object.stokes = stokes; % We take the new stokes no matter what, as we have storage for all 4.
            
            if (isempty(object.radecstokesbeam))
                object.radecstokesbeam = zeros (object.nelem/4, length (object.Theta), length(object.Phi), object.nfreq, length(stokes)); 
            end;
            cohmat = zeros (2,2);
            stI = 0.5*eye(2);                 % 0.5(v_xv_x* + v_yv_y*)
            stQ = stI; stQ (2,2) = -stQ (2,2);% 0.5(v_xv_x* - v_yv_y*)
            stU = rot90(stI);                 % 0.5(v_xv_y* + v_yv_x*)            
            stV = stU; stV(2,2) = -stU(2,1);  % 0.5(v_xv_y* - v_yv_x*)
            stokesmap = containers.Map ({'I','Q','U','V'}, {stI, stQ, stU, stV});
            
            for NF = 1:object.nfreq
                fprintf (1,'%.2f...', object.Freq(NF)/1e6); 
                for NT=1:length(object.Theta)
                    for NP=1:length(object.Phi)
                        for ant = 1:object.nelem/4
                            % tmp = squeeze(object.gain(NT,NP,NF,ant,:,:)); % 2x2 Jones matrix
                            tmp = squeeze (object.vout(ant:ant+1,:,NT,NP,NF));
                            for st = 1:length(stokes)
                                cohmat = tmp*squeeze(stokesmap(upper(stokes(st))))*tmp';
                                % below is already a power due to adding up xx and yy
                                % of the CM.
                                object.radecstokesbeam (ant,NT, NP, NF, st) = cohmat(1,1) + cohmat(2,2);
                            end; 
                         end;
                    end;
                end;
            end;
            fprintf (1, '<-- Generated stokes power beam.\n');           
        end;
        
        % Average RADEC beams averaged over station dipoles
        % Arguments
        %  ignoreants: Antennas to ignore in the averaging
        % Returns:
        %  None: output internally stored
        function  genAvgStationRaDecBeams(object, ignoreants)
            assert (object.simloaded == 1);
            assert (isempty (object.radecstokesbeam) == 0)
            if (isempty (object.radecstnavgbeam) == 1)
                object.radecstnavgbeam = zeros (length(object.station_names), length(object.Phi), length(object.Theta), length(object.Freq), length(object.stokes));
            end;
            tmp = reshape(object.radecstokesbeam, [length(object.station_names),48,length(object.Phi),length(object.Theta),length(object.Freq), length(object.stokes)]);
            for ind = 1:length(stations)
                rem_ants = setdiff(1:48, ignoreants{ind});
                object.radecstnavgbeam(ind,:,:,:,:) = mean (tmp(ind, rem_ants,:,:,:,:), 1);
            end;                        
        end;
            
        % Function to convert RA/DEC beams to l,m. User can specify the
        % frequencies at which to sample the beams, but gets all stokes
        % specified to genRaDecStokesBeam.
        % Arguments:
        %  stokes: Cell array of stokes beams to convert
        %  l,m   : l,m range at which to convert
        %   freq : Array of frequencies at which to generate the power
        %   pattern.
        % Returns:
        %  None: Generated outputs are internally stored.
        function genLMBeams (object, l, m, freq)
            assert (object.simloaded == 1);
            assert (isempty (object.radecstokesbeam) == 0);
            
            if (isempty(l)==1 || isempty(m)==1)
                l = -1:0.01:1;
                m = l;
            end;
            
            assert (min(l) >= -1 & max(l) <= 1);
            assert (min(m) >= -1 & max(m) <= 1);
            if (isempty (freq) == 1)
                % freq = object.Freq(1:3:length(object.Freq)); % Complains when using end;
                freq = object.Freq;
            end;
            
            object.l = single(l); object.m = single(m); % Store user specified l/m positions where the beam is sampled.
            object.usrfreq = freq;
            if (isempty (object.lmstokesbeam) == 1)
                object.lmstokesbeam = single (zeros (object.nelem/4, length(l), length(m), length(freq), length(object.stokes)));
            else
                if (size (object.lmstokesbeam, 2) ~= length (l) || size (object.lmstokesbeam, 4) ~= length(freq))
                    clear object.lmstokesbeam;
                    object.lmstokesbeam = single (zeros (object.nelem/4, length(l), length(m), length(freq), length(object.stokes)));                    
                end;
            end;
            
            [lgrid, mgrid] = meshgrid(l, m);
            dist = sqrt(lgrid.^2 + mgrid.^2);
            thetai = asin(dist);
            thetai(dist >= 1) = 0;
            thetai = thetai * 180/pi; % Convert to degrees.
            phii = mod(atan2(lgrid, mgrid), 2 * pi)*180/pi; % Convert to degrees.

            % Interpolate the simulated beam to the frequency and l,m positions needed 
            % by the user, for all 576 antenna            
            for idx = 1:length(freq)
                for dip = 1:object.nelem/4
                    for st = 1:length(object.stokes)
                        object.lmstokesbeam(dip,:, :, idx, st) = single (interp3 (object.Phi, object.Theta, object.Freq, squeeze(object.radecstokesbeam(dip,:,:,:,st)), phii, thetai, freq(idx) * ones(size(phii))) .* (dist < 1));
                    end;
                end;
            end
            fprintf (1, '<-- Generated interpolated l,m beam.\n');
        end;
       
        % Average LM beams over stations
        % Arguments
        %  ignoreant: Antennas to ignore in the averaging. Specified as a 
        %      cell array of vectors with antenna indices within a station.
        %      Used to eliminate dipoles with the most damaging mutual coupling.
        % Returns:
        %  None: output internally stored
        function  genAvgStationLMBeams(object, ignoreants)
            assert (object.simloaded == 1);
            assert (isempty (object.lmstokesbeam) == 0)
            if (isempty (object.lmstnavgbeam) == 1)
                object.lmstnavgbeam = zeros (length(object.station_names), length(object.l), length(object.m), length(object.usrfreq), length(object.stokes));
            end;
            tmp = reshape(object.lmstokesbeam, [length(object.station_names),48,length(object.l),length(object.m),length(object.stokes), length(object.usrfreq)]);
            for ind = 1:length(object.station_names)
                rem_ants = setdiff(1:48, ignoreants{ind});
                object.lmstnavgbeam(ind,:,:,:,:) = squeeze(mean (tmp(ind, rem_ants,:,:,:,:), 2));
            end;
            
        end;
             
        % Generate the beam for the full AARTFAAC array with choice of
        % stations.
        % Arguments:
        %   stations: Vector of stations to include in the average beam
        % Returns:
        %   None, output is internally stored.
        function genAARTFAACLMBeam (object, stn)
            assert (isempty (object.lmstnavgbeam) == 0);
            if (isempty (stn) == 1)
                stn = 1:6;
            end;
            if (isempty(object.lmaartfaacbeam) == 1)
                object.lmaartfaacbeam = zeros (length (object.l), length (object.m), length(object.usrfreq), length (object.stokes));
            end;
            
            object.lmaartfaacbeam = squeeze(mean (object.lmstnavgbeam(stn,:,:,:,:), 1));
        end;

        % Save the generated beams as .mat files, for e.g., beam model generation.
        % Arguments:
        %    fname : Name of output file. Generated if not provided.
        %    ftype : Type of file. Allowed options:
        %             'mat' : Matlab .mat file (Default)
        %             'fits': FITS file via fitswrite()
        %             'hdf5': HDF5 file via H5write()
        % Returns :
        %    None.
        function saveAARTFAACLMbeam (object, fname, ftype )
            assert (isempty (object.lmaartfaacbeam) == 0);
            if (isempty (ftype) == 1)
                fprintf (2,'saveAARTFAACLMbeam: File type not found, defaulting to .mat');
                ftype = 'mat';
            end;
            if (isempty (fname) == 1)
                fname = sprintf ('%s%s_AARTFAAC_beamshape_%s.%s', object.fsaveprefix, object.arrayconfig, datestr(now(), 30), ftype);
            end;
            switch upper (ftype)
                case 'MAT'
                    save (fname, object.lmaartfaacbeam, object.usrfreq, object.l, object.m);
                case 'HDF5'
                    h5create (fname, '/lmbeamintensity_norm', size (object.lmaartfaacbeam));
                    h5create (fname, '/freq_hz', size (object.usrfreq));
                    h5create (fname, '/l', size (object.l));
                    h5create (fname, '/m', size (object.m));
                    h5info (fname);
                    h5write (fname, '/lmbeamintensity_norm', object.lmaartfaacbeam);
                    h5write (fname, '/freq_hz', object.usrfreq);
                    h5write (fname, '/l', object.l);
                    h5write (fname, '/m', object.m);
                    
                otherwise
                    fprintf (2, 'saveAARTFAACLMbeam: Unknown extension %s. Saving as .mat',ftype);
                    return ;
                
            end;
        end;
        
%%%%%%%%%%% Plotting Related %%%%%%%%%%%%%
        % Setup the plotting surface as necessary
        function fhdl = setupPlot (object, fhdl)
            if (isempty (fhdl))
                object.fhdl = figure ();
            end;
           set (object.fhdl, 'Position', [0, 0, 900, 1100]);
           fhdl = object.fhdl;
        end;
        
        
        % Plot the locations of all dipoles
        % Arguments: 
        %   fhdl: figure handle onto which to plot
        % Returns:
        %   None
        function fhdl = plotPositions (object, fhdl)
            if (isempty(fhdl))
                fhdl = object.setupPlot ([]);
            end;
            assert (object.simloaded == 1);
            
            plot (object.Positions(1,object.dipsel), object.Positions(2,object.dipsel), '.b');
            xlabel ('m'); ylabel ('m');
            title (sprintf ('Dipole positions: %s', object.arrayconfig));
        end;
        
        % Plot a single dipole in a different color among all other dipoles
        % Arguments:
        %   dip: Dipole number to plot
        %  fhdl: Figure handle into which to plot
        function fhdl = plotDipOnPositions (object, dip, fhdl)
            assert (object.simloaded == 1);
            assert (dip < length (object.Positions));
            
            plotPositions (fhdl);
            hold on;
            plot (object.Positions (1,dip), object.Positions (2,dip), 'r*');
            hold off;
        end;
        
        % Image the beampattern of a single dipole
        % Arguments:
        %  fhdl : figure handle to plot to
        %  ant  : Antenna number whose beam needs to be plotted.
        % freqid: Frequency id (in user freq) at which to plot.
        function fhdl = imgEmbedBeam (object, fhdl, ant, freqid)
            if (isempty (fhdl) == 1)
                  figure (fhdl);
            end;      
            assert (isempty (object.l) == 0)
            if (isempty(freqid))
                freqid = 1;
            end;
            norm_fact = max(max(object.lmstokesbeam(ant,:,:,freqid, 1)));
            imagesc (object.l, object.l, object.lmstokesbeam(ant,:,:,freqid,1)./norm_fact);
            colorbar;
            set (gca, 'YDir', 'Normal'); % To match orientation with station images
            set (gca, 'XDir', 'Reverse'); % To match orientation with station images
            text (-0.70, -0.75, sprintf ('%s', object.station_names{station}), 'Color', [1,1,1]);
            title (sprintf ('%s: %f MHz', object.arrayconfig, object.usrfreq(freqid)));
        end;
        
        % Plot the mean primary beam per station over frequency
        % Arguments:
        %   fhdl : Figure handle on which to plot
        %   fnameprefix: Filename prefix with which to store the generated
        %   image.       
        function fhdl = imgStnMeanBeam (object, fhdl, fnameprefix)
            if (isempty(fhdl) == 1)
                fhdl = figure();
            end;
            figure(fhdl);
            
            if (isempty(fnameprefix) == 1)
                fnameprefix = 'station_stI';
            end;
            
            if (isempty (object.lmstnavgbeam) == 1)
                object.genAvgStationLMBeams({[],[],[],[],[],[]});
            end;
            
            for idx = 1:length(object.usrfreq)       
                for station = 1:length(object.station_names);               
                    % subplot (2,3,2*station-1);
                    subplot (3,2,station);
                    b1 = abs(squeeze(object.lmstnavgbeam(station, :,:,idx, 1)));
                    b1 = b1 ./ max(max(abs(b1)));
                    imagesc (object.l,object.m, b1); colorbar;

                    c = contourc (object.l, object.m, b1, [0.75, 0.5, 0.25]);
                    hold on;
                    c1 = c(2,1); c2 = c(2,c1+2);
                    plot (c(1,2:c1+1), c(2,2:c1+1), '-r');
                    plot (c(1, c1+3:c1+c2+2), c(2,c1+3:c1+c2+2), '-k');
                    plot (c(1, c1+c2+4:end-2), c(2,c1+c2+4:end-2), '-c');
                    %axis equal
                    % axis tight
                    set (gca, 'YDir', 'Normal'); % To match orientation with station images
                    set (gca, 'XDir', 'Reverse'); % To match orientation with station images
                    text (-0.67, -0.75, sprintf ('%s', object.station_names{station}), 'Color', [1,1,1]);            
                    hold off;                   
                end;               
                mtit (sprintf ('Stokes-I power patterns, %0.2f MHz, %s', object.usrfreq(idx)/1e6, object.arrayconfig), 'xoff', -0.1, 'yoff', 0.025);
                print ('-depsc', '-r200', sprintf ('%s_%s_%2.0fMHz.eps', fnameprefix, object.arrayconfig, object.usrfreq(idx)/1e6));
                % saveas (gcf, sprintf ('sterp_stations_stokesI_%0.2fMHz.png', object.Freq(idx)/1e6));
                pause;     
                clf;
            end;        
            fprintf (1, '<-- Plotted mean primary beam per station.\n');
        end;
        
        % Plot the aartfaac full beam response over frequency
        % Arguments:
        %   None
        function fhdl = imgAARTFAACBeam(object, fhdl, fnameprefix)
            if (isempty (fhdl) == 1)
                fhdl = figure();
            end;
            assert (isempty (object.lmaartfaacbeam) == 0);
            
            if (isempty(fnameprefix) == 1)
                fnameprefix = 'AARTFAAC_stI';
            end;
            
            for ind = 1:9:length (object.usrfreq)
                for pl = 1:9
                    if ((ind-1)+pl > length(object.usrfreq))
                        break;
                    end;
                    subplot (3,3,pl);
                    norm_fact = max(max(real(object.lmaartfaacbeam(:,:,(ind-1)+pl))));
                    b1 = real(object.lmaartfaacbeam (:,:,(ind-1)+pl))./norm_fact;
                    imagesc (object.l, object.m, squeeze (b1));
                    c = contourc (object.l, object.m, b1, [0.5, 0.25]);
                    hold on;
                    c1 = c(2,1); c2 = c(2,c1+2);
                    plot (c(1,2:c1+1), c(2,2:c1+1), '-r');
                    plot (c(1, c1+3:c1+c2+2), c(2,c1+3:c1+c2+2), '-k');
                    
                    set (gca, 'YDir', 'Normal'); % To match orientation with station images
                    set (gca, 'XDir', 'Reverse'); % To match orientation with station images
                    colorbar;
                    title (sprintf ('%.2f MHz', object.usrfreq((ind-1)+pl)/1e6));
                end;
                mtit (sprintf ('%s Stokes-I average beam over %s.\n', object.arrayconfig, char(object.station_names))); 
                print ('-depsc', '-r200', sprintf ('%s_%s_%2.0fMHz.eps', fnameprefix, object.arrayconfig, object.usrfreq(ind)/1e6));                
                input ('Enter a key to move to the next set of plots');
                clf;
            end;
        end;
        
        % Image the spectral response of the beam per station,
        % unnormalized.
        % Arguments: 
        %  
        % Returns:
        %   Nonezeros (object.nelem/4, length(l), length(m), length(freq), length(stokes));
        
        function imgSpectralEleResp(object, fhdl, fnameprefix)
            if (isempty(fnameprefix) == 1)
                fnameprefix = 'AARTFAAC_spectral_ele';
            end;
            for ind = 1:length(object.station_names)
                subplot (3, 2, ind);
                plot (object.usrfreq /1e6, squeeze(object.lmstokesbeam(ind, :, :)), '.-')
                xlabel ('Freq(MHz)'); ylabel ('Normalized linear gain');
                title (sprintf('%s', object.arrayconfig, object.station_names{ind}));
            end;
            mtit (sprintf ('Beam pattern spectral response: %s FoV', num2str(sample_beam)));             
            print ('-depsc', '-r200', sprintf ('%s_%s_%2.0fMHz.eps', fnameprefix, object.arrayconfig, object.usrfreq(ind)/1e6)); 
        end;
        
        
        % Generate the spectral response per station at the zenith and various elevations.
        % Arguments:
        %  ele: The elevations at which to sample the beam, as a fraction
        %  of the Field of view.
        
%         function plotSpectralEle (object, ele)            
%             % The l/m coordinates of the supplied elevation,
%             % converted to pixel coordinates.
%             if (isempty (ele) == 1)
%                 ele = [0.01, 0.5, 0.7, 0.9];
%             end;
%             sample_pix = int32 (length(object.l)/2) + int32(ele*length(object.l)/2);
% 
%             beam_freq = zeros (6, length (object.usrfreq), length(sample_pix));            
%             for station = 1:length(object.station_names)
%                 for ind = 1:length(object.usrfreq)
%                     beam_freq (station, ind, :) = diag(squeeze(abs(object.lmstnavgbeam(station, sample_pix, sample_pix, ind))));                    
%                 end;
%             end;
% 
%             clf;
%             for ind = 1:length(station_names)
%                 subplot (3, 2, ind);
%                 plot (object.usrfreq/1e6, squeeze(beam_freq(ind, :, :)), '.-')
%                 xlabel ('Freq(MHz)'); ylabel ('Normalized linear gain');
%                 title (sprintf('%s, %s', object.station_names{ind}, object.arrayconfig));                
%             end;
%             mtit (sprintf ('Beam pattern spectral response: %s FoV', num2str(sample_beam))); 
%             % saveas (gcf, 'beam_spectral_arrayconfig_response.png');
%         end
        
        % Generate a movie of the generated embedded element responses,
        % with the movie axis over elements for a fixed frequency.
        % Arguments: 
        %  lmbeam  : The lm beam to be plotted.
        %   fhdl   : figure handle onto which to plot.
        % savefname: prefix with which to save plots. Default is the
        %            timestamp.
        %  freqind : Frequency at which to generate the movie.
        %  freq    : Frequencies at which the lm beam is available.
        %    stokes: The stokes for which the beams have to be generated.
        % Returns:
        %   None
        
        function imgMovieEmbeddedResp(object, fhdl, savefname, freqind, stokes)
            % Generate a movie of the stokes-I beam for each antenna's polarized response within a station, for the specified frequency.            
            if (isempty(freqind) == 1)
                freqind = 1;
            end;
                
            assert (freqind <= max(object.usrfreq));
            if (isempty(savefname) == 1)
                savefname = strcat (sprintf ('%s_%s_%s_beam_pattern_%fMHz', datestr(now, 30), object.arrayconfig, object.usrfreq(freqind)), '.avi');
            end;
            vidObj = VideoWriter (savefname, 'Motion JPEG AVI');
            vidObj.FrameRate = 5;
            open (vidObj);
            figure (fhdl);
            for station = length (object.stations)
                for  ind = 1:48 
                    clf
                    ant = (station-1)*48 + ind;
                    stat = (station-1)*48+1;
                    subplot (1,3,1);
                    object.plotPositions(fhdl);
                    object.plotDipkOnPosition (2*ant, hdl);
                    
                    subplot (1,3,2);
                    tmp = squeeze(abs(lmbeam(ant,:,:,freqind)));
                    imagesc (l,l,tmp./max(max(tmp))); colorbar;
                    title (sprintf ('CS00%d %s:Ant %d, %.2f MHz', station+1, object.arrayconfig, ind, freq(freqind)/1e6));

                    subplot (1,3,3); 
                    tmp = squeeze(mean(squeeze(abs(lmbeam(stat:ant,:,:,freqind))), 1)); 
                    imagesc (l,l,tmp./max(max(tmp))); colorbar;
                    title (sprintf ('CS00%d %s Mean beam, %.2f MHz', station+1, object.arrayconfig, freq(freqind)/1e6));
                    pause(0.1);
                    currFrame = getframe (fig1);
                    writeVideo (vidObj, currFrame);
                end;
            end;
            close(vidObj);
        end;
        
        % plot the spectral variation of the full AARTFAAC beam power contours
        % Arguments:
        %  pwrlev: vector of gain levels at which to plot the spectral
        %          response
        function pltSpectralPowerResponse(object, fhdl, pwrlev, fnameprefix)
            col = {'b', 'g',' r','c', 'm', 'y', 'k', 'w'};
            leginfo = {};
            
            if (isempty(pwrlev))
                pwrlev = [0.75, 0.5, 0.25];
            end;
            
            if (isempty(fnameprefix) == 1)
                fnameprefix = 'AARTFAAC_spectral_resp';
            end;
            
            for pind = 1:length(pwrlev)
                clf;
                hold on;
                for ind = 1: length(object.usrfreq)/2
                    b1 = abs(object.lmaartfaacbeam (:,:,2*ind))./max(max(abs(object.lmaartfaacbeam(:,:,2*ind))));
                    c = contourc(object.l, object.l, b1, [pwrlev(pind), pwrlev(pind)]);
                    plot (c(1,2:end), c(2,2:end), char(col(ind)));
%                     c1 = c(2,1); c2 = c(2,c1+2);
%                     plot (c(1,2:c1+1), c(2,2:c1+1), '-r');
%                     plot (c(1, c1+3:c1+c2+2), c(2,c1+3:c1+c2+2), '-k');
%                     plot (c(1, c1+c2+4:end-2), c(2,c1+c2+4:end-2), '-c');
%                     hold off;            
                    leginfo{ind} = [num2str(object.usrfreq(ind)/1e6) 'MHz'];           
                 end;
                 legend (leginfo);
                 title (sprintf ('%s, Gain level %f.\n', object.arrayconfig, pwrlev(pind)));
                 grid on;                                                  
                 hold off;
                 print ('-depsc', '-r200', sprintf ('%s_%s_%2.0flev.eps', fnameprefix, object.arrayconfig, pwrlev(pind)*100));
                 input ('Enter a key to move to next power level');
            end;
        end;

            
        % Phase up all the element embedded beam pattern in the direction dir
        % Arguments:
        %  azi/ele: Direction in the sky at which to point the tile
        %  beam, rad.
        % Returns: 
        %  tilebeam : The average beam over the full station
        % TODO
        function tilebeam = phaseTile (object, azi, el)
            assert (object.lba == 0)
            
            % Extract out the position of the 16 dipoles making up a tile
            % wrt. center dipole.
            % Find the baseline vector between the center of the tile and
            % every other antenna
            uloc = meshgrid (object.Positions(1,1:16)) - meshgrid(object.Positions(1,1:16)).';
            vloc = meshgrid (object.Positions(2,1:16)) - meshgrid(object.Positions(2,1:16)).';
                                   
            % We carry out the weight determination for a random single
            % tile, as all tiles have the same configuration.
            
            % Represent the pointing direction as a vector
            point_x = cos(ele)*cos(phi);
            point_y = cos(ele)*sin(phi);
            
            % Baseline phase is the dot product of the baseline vector and
            % the pointing vector.
            bline_ph = dot ([uloc(:) vloc(:)], [point_x, point_y]);
            bline_amp =hypot (uloc(:), vloc(:));
            
            % Create a weight vector with the directional phase
            % information. Since this is identical for every tile, it can
            % also be applied directly to signals from all tiles.
            
            
            

            % Carry out beamforming per tile in the pointing direction,.
            % independently for each tile 
        end;
        
       

%%%%%%%%%%%%%% Tasks %%%%%%%%%%%%%
        % Function to create an average AARTFAAC beam 
        % Arguments:
        %   freq: Frequencies at which to generate the beam pattern
        % Returns:
        %  None
        function createAvgbeam (object, stokes, l, m, freq)
            fprintf (2, '# Stokes currently ignored, only stokes-I beam is generated.\n');
            assert (object.simloaded == 1);
            if (isempty(freq) == 1)
                freq = object.Freq; % object.Freq(1:3:length(object.Freq));
            end;
            assert (min(freq) >= min(object.Freq));
            assert (max(freq) <= max(object.Freq));
            fprintf (1, '<-- Generating output voltage into vout member...\n');
            object.genOutputVoltage();
            
            % fprintf (1, '<-- Generating Jones matrices into gains member...\n');
            % object.genJonesMat ();
            
            fprintf (1, '<-- Generating Stokes beams in RA/Dec into radecstokesbeam member...\n');
            object.genRaDecStokesBeam([]);
            
            fprintf (1, '<-- Generating Stokes embedded beams in l,m into lmstokesbeam member...\n');
            object.genLMBeams (l, m, freq);

            fprintf (1, '<-- Generating station average Stokes beams in l,m into lmstnavgbeam member...\n');
            object.genAvgStationLMBeams ({[],[],[],[],[],[]});
            
            fprintf (1, '<-- Generating AARTFAAC average Stokes beams in l,m into lmaartfaacbeam member...\n');
            object.genAARTFAACLMBeam ([]);
            
        end;
        
        function genPlots (object)
            % Setup plot window
            fhdl = object.setupPlot(object.fhdl);

            % -- Beam patterns averaged over all 6 AARTFAAC stations, as a func. of
            % freq.
            object.imgAARTFAACBeam(object.fhdl, []);

            % -- Per station average stokes-I power patterns over frequency
            object.imgStnMeanBeam (object.fhdl, []);

            % -- Spectral response (contour plot) of average stokes-I power pattern at 
            %    fixed gains.
            object.pltSpectralPowerResponse (object.fhdl, [], []);
        end;
            

    end; % End of methods
    
end % End of classdef

