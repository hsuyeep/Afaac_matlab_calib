% Script to simulate visibilities from a specific array configuration, and source locations.
% pep/28Apr15

array_rad = 5; % meters
Nelem2rad = 100;
array_spacing = array_rad/Nelem2rad;
arraysampling = [-array_rad:array_spacing:array_rad]; % Array extent = -5:5m, elements of array placed every 0.5m (=lambda/2 for 1m lambda).

% Create a rectangular grid of antenna, equi-spaced
[xpos, ypos] = meshgrid (arraysampling);

% Create a disk of antennas, with radius being the extent of array sampling.
discsel = (sqrt(xpos.^2 + ypos.^2) < array_rad);
xpos_disc = xpos (discsel);
ypos_disc = ypos(discsel);

% Create a ring of antenna, with inner rad = inner_rad
inner_rad = 3; % meters
innerdiscsel = (sqrt(xpos.^2 + ypos.^2) < inner_rad);
ringsel = discsel - innerdiscsel;
xpos_ring = xpos(ringsel==1);
ypos_ring = ypos(ringsel==1);


Nelem = length (xpos(:));

% Show the array layout. A different colored dot for each row of elements in the array.
plot (xpos, ypos, '.'); 

l0 = 0, m0 = 0; % Locate sources of unit amplitude at the locations (in l,m units) as specified above.

we = exp (2*pi*i*(xpos(:)*l0 + ypos(:)*m0)); % Generate phasor due to location of each element
V = we * we'; % Generate the visibilities for the system at hand.

img_l = [-1:0.01:1];
img_m = img_l;

% Create an image using acm2skymap:
map = acm2skyimage (V, xpos(:), ypos(:), 299792458, img_l, img_m);

imagesc (img_l, img_m, real(map)); colorbar; axis tight;
xlabel ('l'); ylabel ('m');

