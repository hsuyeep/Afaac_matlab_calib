% Function to carry out imaging of the CCM using an FFT.
% ccm    = cross correlation matrix of observations
% xpos   = x-positions of antennas in ITRF frame
% ypos   = y-positions of antennas in ITRF frame
%         The (u,v) values are generated from xpos and ypos by taking their
%         differences, which give baseines
% uvmax  = The max. uv range (in m) of the generated grid. To form the
% hermitiean symmetric uv grid, the range is -uvmax:uvmax

% gridsize = resolution of UV grid (in m).
% visgrid = gridded UV values after interpolation of observed visibilities
% map    = fourier inverted visibilities
% weights = size (u) x size (v) matrix containing a selection matrix
function [map, weights, visgrid] = fft_imager (ccm, xpos, ypos, uvmax, gridsize)

%% create the observed uv values
uobs = meshgrid (xpos) - meshgrid (xpos).';
uobsvec = uobs (:);
vobs = meshgrid (ypos) - meshgrid (ypos).';
vobsvec = vobs (:);
% umax = max (uobs (:)); umin = min (uobs (:));
% vmax = max (vobs (:)); vmin = min (vobs (:));

% create a grid with specified settings
ugrid = [-uvmax:gridsize:uvmax];
npoints = length (ugrid);

% Assuming u and v spacing is identical
visgrid = zeros (length (ugrid));
weights = visgrid;

%% Grid the observed visibilities.
for ind = 1:length (ccm)
    % for every observed visibility, decide on which grid point it should
    % fall
    uobsind  = uobsvec(ind)/gridsize + npoints/2; % uvmax;
    vobsind  = vobsvec(ind)/gridsize + npoints/2; % uvmax;
    % ugridind = fix (uobs(ind)/gridsize) + uvmax ; % Make  indices positive
    % vgridind = fix (vobs(ind)/gridsize) + uvmax ; % for Matlab
    uvamp    = abs (ccm(ind));
    uvphasor = ccm(ind)/uvamp;
    
    % For distributing the amplitude of the observed visibility such that
    % the center of mass lies at the observed position
    % First find the center of mass along the u axis, but at observed v
    % position
    ugridh = ceil (uobsind);
    ugridl = floor (uobsind);
    %sig_uh = ((uobsind - ugridh) * uvamp)/gridsize;
    %sig_ul = uvamp - sig_uh; 
    sig_uh = abs ((uobsind - ugridh) * uvamp); %/gridsize;
    sig_ul = abs ((uobsind - ugridl) * uvamp); % uvamp - sig_uh; 
    
    
    % Then distribute the u axis low and high ampl. along the v axis 
    vgridh = ceil (vobsind);
    vgridl = floor (vobsind);
    %sig_uh_vh = ((vobsind - vgridh) * sig_uh)/gridsize;
    %sig_uh_vl = sig_uh - sig_uh_vh;
    %sig_ul_vh = ((vobsind - vgridl) * sig_ul)/gridsize;
    %sig_ul_vl = sig_ul - sig_ul_vh;
    
    sig_uh_vh = abs ((vobsind - vgridh) * sig_uh); %/gridsize;
    sig_uh_vl = abs ((vobsind - vgridl) * sig_uh); % sig_uh - sig_uh_vh;
    sig_ul_vh = abs ((vobsind - vgridh) * sig_ul); %/gridsize;
    sig_ul_vl = abs ((vobsind - vgridl) * sig_ul); % - sig_ul_vh;
    
    % Generate interpolated visibility for this grid point
    visgrid (ugridh, vgridl) = visgrid (ugridh, vgridl) + sig_uh_vl * uvphasor;
    visgrid (ugridh, vgridh) = visgrid (ugridh, vgridh) + sig_uh_vh * uvphasor;
    visgrid (ugridl, vgridl) = visgrid (ugridl, vgridl) + sig_ul_vl * uvphasor;
    visgrid (ugridl, vgridh) = visgrid (ugridl, vgridh) + sig_ul_vh * uvphasor;
    
    weights (ugridh, vgridl) = weights (ugridh, vgridl) + 1;
    weights (ugridh, vgridh) = weights (ugridh, vgridh) + 1;
    weights (ugridl, vgridl) = weights (ugridl, vgridl) + 1;
    weights (ugridl, vgridh) = weights (ugridl, vgridh) + 1;
end

%%
map = fft2 (visgrid);
%end