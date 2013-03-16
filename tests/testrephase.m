% Script to test out imaging with rephasing.
% Idea is to generate a sequence of rephased acms, and their images, to 
% understand if rephasing is occuring correctly or not.
% NOTE: Uncomment section with 'keyboard' to enable stepping through
% the mosaic creation mechanism.
% pep/02Nov12

close all; clear all;

% Load antenna positions
load ('poslocal.mat');
load ('srclist3CR.mat');
uloc = meshgrid (poslocal(:,1)) - meshgrid (poslocal(:,1)).';
vloc = meshgrid (poslocal(:,2)) - meshgrid (poslocal(:,2)).';

% imaging parameters
nfacet = 3;      % Create nfacet x nfacet facet images across FoV;
facetsize = 256; % size of the facet in pixels 
Nuv = 700;
duv = 2.5; uvpad = facetsize *nfacet; % Each facet is 512x512 pixels

%% load an acm
% fid = fopen ('~/WORK/AARTFAAC/rawdat/SB000_ch30-35_5sec_3hr_cal.bin', 'rb');
fid = fopen ('~/WORK/AARTFAAC/Reobs/11Jul12/LBA_OUTER_BAND60/SB000_LBA_OUTER_1ch_cal.bin', 'rb');
[acc, tobs, freq] = readms2float (fid, -1, -1);
if (isempty(tobs) == true)
	disp ('End of file reached!'); 
	return;
end;

% rephase!
% 2 units = full sky => 1 unit shift in l/m is 90deg, or 320 mins in time.
% If we want to generate sky images covering, say, 80% of the FWHM of 120deg,
% then we cover 96deg, or ~320 mins.
m_reph = 0; l_reph = 0;

msize = nfacet*facetsize; % One side of the mosaic.

lmtopix = 2/uvpad; % size in l,m units of each pixel in generated mosaic.
moff = msize / nfacet; % Offset of center of facet.

mosaic = zeros (msize);
prevmosaic = mosaic;
laxis = zeros (1, nfacet*facetsize);
zenimg = figure;
mosimg = figure;
diffimg = figure;
pro = figure;
mask = [];

%% Start generating mosaic images for different timeslices.
while (isempty(acc) == false)
lind = 2; 
 for l_reph = -1+1/nfacet:2/nfacet:1-1/nfacet
 mind = 2;
 % NOTE: Unconventional indexing on lind needed to match orientation
 % of image wrt. original.
 laxis ((lind-2)*(-moff)+1:(lind-3)*(-moff)) = [l_reph-(facetsize/2-1)*lmtopix:lmtopix:l_reph+facetsize/2*lmtopix];
 
 for m_reph = -1+1/nfacet:2/nfacet:1-1/nfacet
  disp (['l/m_reph: ' num2str([l_reph, m_reph])]);
  if ((l_reph*l_reph + m_reph*m_reph) < 1)
    [reph_acc] = rephaselm (acc, l_reph, m_reph, tobs, freq, posITRF);
    [radecmap1, calmap, calvis, l, m] = fft_imager_sjw_radec (reph_acc (:), ... 
                    uloc(:), vloc(:), duv, Nuv, uvpad, tobs, freq, 0);

    if isempty (mask)
     mask = NaN (length (l));
     mask ((meshgrid (l).^2 + (meshgrid(l).').^2) < 1) = 1;
    end;
    
    map = abs (calmap) .* mask;
    [zen_dr, zen_sig] = getimagedr (map, 64, 3);
    imgmax = max(max(map));
    imgmin = min(min(map));     
    if (lind == 1  & mind == 1) % Comment the conditional to see every rephased image.
%      figure (zenimg);
%      imagesc (l, m, map);
%      colorbar;
%      caxis ([imgmin 0.5*imgmax]);   
%      title (sprintf ('Time: %f,Phase center: l=%f, m=%f. DR = %4.1f, sig = %4.1f', ... 
%             tobs, l_reph, m_reph, zen_dr, zen_sig));
%      xlabel ('East - West');
%      ylabel ('North - South');
      zenithmap = map;
    end
    mosaic (lind*moff+1:(lind+1)*moff, mind*moff+1:(mind+1)*moff) = ...
         map (uvpad/2-facetsize/2+1:uvpad/2+facetsize/2, uvpad/2-facetsize/2+1:uvpad/2+facetsize/2);
    disp (['(l,m) range: (' num2str([l_reph-facetsize*lmtopix,l_reph+facetsize*lmtopix]) ... 
'), (' num2str([m_reph-facetsize*lmtopix,m_reph+facetsize*lmtopix]) ')']);
    
    
% Uncomment to see mosaic being created with every facet.
%     figure (mosimg);
%     imagesc (mosaic);
%     pause;

    mind = mind - 1;
  end
 end
 lind = lind - 1;
end

% Create an appropriate mask for the mosaic image.
mosmask  = zeros (length(laxis));
mosmask ((meshgrid(laxis).^2 + (meshgrid(laxis).').^2) < 1) = 1;
mosaic = mosaic .* mosmask;

% display the mosaic.
figure (mosimg);
[mos_dr, mos_sig] = getimagedr (mosaic, 64, 3);
imagesc (laxis, laxis, mosaic);
mosmax = max(max(mosaic(mosmask == 1)));
mosmin = min(min(mosaic(mosmask == 1)));
title (sprintf ('Mosaic: %d facets (%d pix). DR = %4.1f, sig = %4.1f (max: %f, min: %f', ...
       nfacet, facetsize, mos_dr, mos_sig, mosmax, mosmin));
colorbar;
caxis ([mosmin 0.5*mosmax]);             
xlabel ('East - West');
ylabel ('North - South');
[sel, sel_l, sel_m] = overplot3cr (tobs, srclist3CR, 150, mosimg);
     

% display the difference between the zenith image and the mosaic
figure (zenimg); % NOTE: Using zenith image plot for this!
zenithmap (mask == 0) = 1;
mosaic (mosmask == 0) = 1;
zen_mos_diff = (abs(mosaic) ./ abs(zenithmap));
mosmax = max(max(zen_mos_diff(mosmask == 1)));
mosmin = min(min(zen_mos_diff(mosmask == 1)));
imagesc (laxis, laxis, zen_mos_diff);
title (sprintf ('Mosaic ratio: max = %f, min = %f', ...
       mosmax, mosmin));
colorbar;
caxis ([0 2]);

% Difference between timeslices
mosdiff = mosaic - prevmosaic;
prevmosaic = mosaic;
diffmean = mean (mosdiff(:));
diffsig = std (mosdiff(:));
diffmax = max(max(mosdiff(mosmask == 1)));
diffmin = min(min(mosdiff(mosmask == 1)));
figure(diffimg)
imagesc (laxis, laxis, mosdiff);
title (sprintf ('Mosaic Time diff: max = %f, min = %f, mean = %f, var = %f', ...
       diffmax, diffmin, diffmean, diffsig));
% caxis ([0 2]);
colorbar;
xlabel ('East - West');
ylabel ('North - South');



% display mosaic and time difference profiles, and their transforms
figure (pro);
l_pro = sum (mosaic, 1);
m_pro = sum (mosaic, 2);
l_prodiff = sum (mosdiff, 1);
m_prodiff = sum (mosdiff, 2);

subplot (2, 5, 1)
plot (laxis, l_pro);

subplot (2, 5, 2)
L_pro = fft(l_pro);
plot (abs(L_pro(2:end/2)));

subplot (2, 5, 3)
plot (laxis, m_pro);

subplot (2, 5, 4)
M_pro = fft(m_pro);
plot (abs(M_pro(2:end/2)));

subplot (2, 5, 5);
hist (mosaic (:), 100);

subplot (2, 5, 6);
plot (laxis, l_prodiff);

subplot (2, 5, 7)
L_prodiff = fft (l_prodiff);
plot (abs(L_prodiff(2:end/2)));

subplot (2, 5, 8)
plot (laxis, m_prodiff);

subplot (2, 5, 9)
M_prodiff = fft(m_prodiff);
plot (abs(M_prodiff(2:end/2)));

subplot (2, 5, 10)
hist (mosdiff(:), 100);

 % Get next timeslice.
 [acc, tobs, freq] = readms2float (fid, -1, -1);
 if (isempty(tobs) == true)
 	disp ('End of file reached!'); 
	return;
 end;
 pause;
end;

fclose (fid);
