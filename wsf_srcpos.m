% Function to estimate model src positions from the ACM.
% pep/05Oct12
% Arguments:
%         acc : Uncalibrated ACM
%        cal1 : Direction independent gains per antenna element.
%	     freq : Frequency of ACM.
%     Sigman1 : Noise estimate from calibration of ACM.
%        nsrc : Number of source for which position needs to be estimated.
%      rodata : Read only data containing constants.
%     phisrc0 : Initial estimates of azimuth angles of all skymodel sources.
%      thsrc0 : Initial estimates of elevation angles of all skymodel sources.
%         sel : selection vector to choose relevant srcs from all skymodel srcs.
%         opt : Optimization parameters, created by caller using 'optimset'
%
% Returns:
%   th/phisrc_cat : Catalog positions of selected sources at given time
%   th/phisrc_wsf : Estimated positions of selected sources.
%   fval/out   : Function evaluation value, output structure from fminsearch.
% ---------  Estimate model source positions using WSF algorithm ----------
% Ref     : Viberg et. al, IEEE Trans. sig. proc, vol. 39, No. 11, 1991

function [thsrc_cat, phisrc_cat, thsrc_wsf, phisrc_wsf, fval, out] = ... 
    wsf_srcpos (acc, cal1, freq, Sigman1, nsrc, rodata, phisrc0, thsrc0, ...
                sel, opt, debug)

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

    % Construct weight matrix which gives lowest asymptotic variance 
	% (see eq. 25 of [Ref].)
    Wopt = (diag(d(1:nsrc)) - mean(diag(squeeze(Sigman1))) * eye(nsrc))^2 / ... 
           diag(d(1:nsrc));

    % Extract signal eigenvectors corresponding to top nsrc eigen values.
    Es = v(:, 1:nsrc);  
    EsWEs = Es * Wopt * Es'; % WSF weight matrix.

	% Carry out maximization of the log likelihood function.
    [theta, fval, exitflag, out] = ... 
 			fminsearch(@(x) WSFcostfunITRF(x, EsWEs, diag(conj(1 ./ cal1)), ...
                       freq, rodata.posITRF_fl), [phisrc0; thsrc0], opt);
    phisrchat = zeros (length(rodata.srcsel),1); 
    phisrchat(sel) = theta(1:nsrc);
    thsrchat  = zeros (length(rodata.srcsel),1); 
    thsrchat (sel) = theta(nsrc+1:end);

    if (debug > 0)
		disp('Position estimation using WSF done: ');
    end 
    if (debug > 1)
		disp ('Catalog positions: '); disp ([thsrc0 phisrc0]');
		disp ('WSF     positions: '); disp ([thsrchat phisrchat]');
    end
    thsrc_cat  = zeros (length(rodata.srcsel), 1); thsrc_cat (sel) = thsrc0;
    phisrc_cat = zeros (length(rodata.srcsel), 1); phisrc_cat (sel) = phisrc0;
    thsrc_wsf  = thsrchat;
    phisrc_wsf = phisrchat;
