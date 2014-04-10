% Program to carry out wide field imaging using facets.
% pep/05Nov12
% Arguments:
% 	acc  : calibrated ACM to image
%	tobs : Time of this observation
%	freq : Freq. of this observation (Hz)
% poslocal: antenna positions in ITRF coords, m, local to CS002
% nfacets: The number of facets which will make up the full field of view.
%		 : in one direction.
% pix2facet: The number of pixels making up each facet.
function [image] = facet(acc, tobs, freq, nfacets, pix2facet, poslocal)
	C = 299792458; % m/s

	% Grid the sky based on the number of facets, and determine facet centers.
	% Full field of view = -2 in direction cosines.
	lmrange = abs(-2/nfacets);  % Extent of each facet.
	l_range = [-1+lmrange/2:lmrange:1-lmrange/2];
	m_range = l_range;
	facet_dl = lmrange/pix2facet; % Resolution of each facet image pixel.
	facet_duv = freq/(C*lmrange);

	uloc = meshgrid (poslocal(:,1)) - meshgrid (poslocal(:,1)).';
	vloc = meshgrid (poslocal(:,2)) - meshgrid (poslocal(:,2)).';

	% Gridding parameters; added on 11Apr14	
	gparm.type = 'pillbox';
	gparm.fft = 1;    % Carry out FFT imaging
    gparm.duv = facet_duv;% Default image just the full Fov (-1<l<1)
    gparm.Nuv = 256;    % size of gridded visibility matrix
    gparm.uvpad = pix2facet;  % specifies if any padding needs to be added

	% For each facet
	for l_ind = 1:nfacets
		for m_ind = 1:nfacets

		% Carry out rephasing of the array towards facet center
		reph_acc = rephaselm (acc, l_range(l_ind), m_range(m_ind), tobs, freq,...
							poslocal);

		disp (['Imaging parameters: ' num2str([facet_duv, pix2facet])]);

		% Carry out recomputing of the uvw coordinates towards the facet center.
		% TODO
	
		% carry out imaging of the facet
  		[radecmap1, calmap, calvis, l, m] =  ...
				fft_imager_sjw_radec (reph_acc (:), ... 
                    uloc(:), vloc(:), gparm, [], [], tobs, freq, 0);

		figure;
		imagesc (abs(calmap));
		end;
	end;
