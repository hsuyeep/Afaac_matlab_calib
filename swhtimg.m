% Script to attempt spherical harmonic imaging from AARTFAAC.
% pep/21May15

function [swht, blm] = swhtimg()
fname = '~/WORK/AARTFAAC/Reobs/20Nov13/r02/SB002_LBA_OUTER_8b2sbr02_1ch.bin';
load ('poslocal_outer.mat', 'poslocal');
doplot = 0;

Nelem = 288;
fid = fopen (fname, 'rb');
[acc, t_obs, freq] = readms2float (fid, 100, -1, Nelem);

% Single station ACM.
acc = acc (1:48, 1:48);
k0 = (2*pi*freq)/299792458; % Wavenumber for this observation.

lmax = 2; % int32(Nelem/50);
swht = zeros (lmax+1, 2*lmax+1);
acc = acc (:);

% Spherical coordinates of each baseline
uloc = meshgrid (poslocal(:,1)) - meshgrid (poslocal(:,1)).';
vloc = meshgrid (poslocal(:,2)) - meshgrid (poslocal(:,2)).';
wloc = meshgrid (poslocal(:,3)) - meshgrid (poslocal(:,3)).';

[blth, blphi, blR] =  cart2sph (uloc(:), vloc(:), wloc(:)); % ITRF m, local to CS002.
blR = blR*freq/299792458; % Radial component in lambda units.

if (doplot == 1)

end;

% Compute the spherical harmonic transform coefficients for each 
for l = 0:lmax
	fprintf (2, '.');
	for m = -l:l
		mabs = abs(m);
		for ind = 1:length(acc)
			% Generate the spherical bessel function for these quantal numbers.
			if (blR(ind) == 0)
				if (l == 0) 
					bes = 1;
				else
					bes = 0;
				end;
			else
				bes = sqrt (pi./(2*blR(ind)*k0)).*besselj(l+0.5, blR(ind)*k0);
			end;

			% Generate the spherical harmonic coefficients for these quantal numbers
			Plm = legendre(l,cos(blth(ind)));
			if l~=0
				Plm = squeeze(Plm(mabs+1,:,:));
			end
			a1 = ((2*l+1)/(4*pi));
			a2 = factorial(l-mabs)/factorial(l+mabs);
			C = sqrt(a1*a2);
			Ylm = C*Plm.*exp(i*mabs*blphi(ind));
			if (m < 0)
				Ylm = conj(Ylm)*(-1)^mabs;	
			end;
			

			swht (l+1,m+l+1) = swht (l+1,m+l+1) + acc(ind)*bes*Ylm;
		end;
	end;
end;

fprintf (1, 'Computed shwt coe.\n');
% Compute the brightess
for ind = 1:size (swht,1)
	blm(ind, :) = swht(ind, :)/(4*pi)*(-i)^(-ind);
end;


% Reproject to local coordinates
res = 0.1;
lloc = meshgrid ([-1:res:1]);
mloc = lloc';
[phisph, rsph] = cart2pol(lloc, mloc);
thsph = asin (rsph);
img = zeros (size (lloc));


for xpix = 1:size (img,1)
	fprintf (2, '+');
	for ypix = 1:size (img,2)
		if (rsph(xpix, ypix) <= 1)
			for l = 0:lmax
				for m = -l:l
					mabs = abs(m);
					% Generate the spherical harmonic coefficients for these quantal numbers
					Plm = legendre(l,cos(thsph(xpix,ypix)));
					if l~=0
						Plm = squeeze(Plm(mabs+1,:,:));
					end
					a1 = ((2*l+1)/(4*pi));
					a2 = factorial(l-mabs)/factorial(l+mabs);
					C = sqrt(a1*a2);
					Ylm = C*Plm.*exp(i*mabs*phisph(xpix,ypix));
					if (m < 0)
						Ylm = conj(Ylm)*(-1)^mabs;	
					end;
			
					img(xpix, ypix) = img(xpix, ypix) + blm(l+1,m+l+1)*Ylm;	
				end;
			end;
		end;
	end;
end;


