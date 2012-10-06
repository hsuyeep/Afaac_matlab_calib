% Function to estimate model src positions from the ACM.
% pep/05Oct12
% Arguments:
%   uloc/vloc : u/v positions of the elements in ITRF local cordinates
%   flagant   : vector containing the antennas to be flagged, as numbers 
%               between 1-Nant
% Returns:
%   u/vloc_flag: vectors containing the u/v positions of the remaining 
%				unflagged elements. 
% ---------  Estimate model source positions using WSF algorithm ----------
% Ref     : Viberg et. al, IEEE Trans. sig. proc, vol. 39, No. 11, 1991

function [thsrc_cat, phisrc_cat, thsrc_wsf, phisrc_wsf] = ... 
    wsf_srcpos (acc, cal1, freq, Sigman1, nsrc, posITRF_fl, phisrc0, thsrc0, ...
                sel, srcsel, debug)

    % whitening of the array covariance matrix (required for DOA estimation)
    % acc = acc ./ sqrt(diag(acc) * diag(acc).');

    % Carry out an Eigenvalue decomposition of ACM.
    [v, d] = eig(acc);

    % Sort eigen values in decreasing order. The vector 'order' holds the index
    % in original 'd' vector of eigenvalues in sorted order. Store result back 
    % into 'd'.
    [d, order] = sort(diag(abs(d)), 'descend');

    % Reorder the eigenvectors as per sorting order.
    v = v(:, order);

    % Construct weight matrix which gives lowest asymptotic variance (see ref.)
    Wopt = (diag(d(1:nsrc)) - mean(diag(squeeze(Sigman1))) * eye(nsrc))^2 / ... 
           diag(d(1:nsrc));

    % Extract signal eigenvectors corresponding to top nsrc eigen values.
    Es = v(:, 1:nsrc);  
    EsWEs = Es * Wopt * Es'; % WSF weight matrix.

	% Carry out maximization of the log likelihood function.
    theta = fminsearch(@(x) WSFcostfunITRF(x, EsWEs, diag(conj(1 ./ cal1)), ...
                       freq, posITRF_fl), [phisrc0; thsrc0]);
    phisrchat = zeros (length(srcsel),1); 
    phisrchat(sel) = theta(1:nsrc);
    thsrchat  = zeros (length(srcsel),1); 
    thsrchat (sel) = theta(nsrc+1:end);

    if (debug > 0)
		disp('Position estimation using WSF done: ');
    end 
    if (debug > 1)
		disp ('Catalog positions: '); disp ([thsrc0 phisrc0]');
		disp ('WSF     positions: '); disp ([thsrchat phisrchat]');
    end
    thsrc_cat  = zeros (length(srcsel), 1); thsrc_cat (sel) = thsrc0;
    phisrc_cat = zeros (length(srcsel), 1); phisrc_cat (sel) = phisrc0;
    thsrc_wsf  = thsrchat;
    phisrc_wsf = phisrchat;

    % srcposhat = [cos(phisrchat) .* cos(thsrchat), ... 
	%              sin(phisrchat) .* cos(thsrchat), ... 
    %              sin(thsrchat)];
