% Script to extract the noise floor from an image.
% Subdivides the image into an MxN grid, extracts a
% separate noise value from each grid section after 
% carrying out a sigma clipping of the sources in each
% grid section.
% pep/22Oct15

% Arguments:
%   img : Image array to operate on
% thresh: Sigma threshold to clip on
%  M,N  : Sections into which to grid the image.
% Returns:
%   rms : Matrix of size M,N, each entry has the robust rms of the 
%         pixels of that region.
%   mu  : Matrix of size M,N, each entry has the robust mean of the 
%         pixels of that region.
% selpix: Matrix of size M,N with the total number of selected pixels used.
% pix2grid: Scalar containing the number of pixels in each grid cell.

function [rms, mu, selpix, pix2grid] = genimgnoisefloor (img, thresh, M, N) 
    assert (isempty(img) == 0);
    if (isempty (thresh)) thresh = 3; end;
    if (isempty (M)) M = 3;  end;% Default is a 3x3 grid.
    if (isempty (N)) N = 3;  end;

    rms = zeros (M,N);
    mu = zeros (M,N);
    selpix = zeros (M,N);
    rstride = floor(size (img, 1)/M);
    cstride = floor(size (img, 2)/N);
    pix2grid = [rstride, cstride];
    for rind = 1:M
        for cind = 1:N
            ro = (rind-1)*rstride + 1;
            co = (cind-1)*cstride + 1;

            % Extract out the valid pixels from each grid region
            pix = img(ro:ro+rstride, co:co+cstride);

            % Carry out an iterative sigma clipping of each grid
            [m,v,sel] = robustmean (pix(:), thresh);

            rms(rind, cind) = v;
            mu(rind, cind) = m;
            selpix(rind,cind) = sum(sel);

        end;
    end;
    % 
