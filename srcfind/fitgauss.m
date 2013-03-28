% Script to fit a gaussian to the provided data.
% Adapted from the TKP codebase.
% pep/20Dec12
%
%    Calculate source positional values by fitting a 2D Gaussian
%    Args:
%
%        data (numpy.ndarray): Pixel values
%
%        params (dict): initial fit parameters (possibly estimated
%            using the moments() function, above)
%
%    Kwargs:
%
%        fixed (dict): parameters & their values to be kept frozen (ie, not
%            fitted)
%
%        maxfev (int): maximum number of calls to the error function
%
%    Returns:
%
%        (dict): peak, total, x barycenter, y barycenter, semimajor,
%            semiminor, theta (radians)
%
%    Raises:
%
%        ValueError (in case of a bad fit)
%
%    Perform a least squares fit to an elliptical Gaussian.
%    If a dict called fixed is passed in, then parameters specified within the
%    dict with the same names as fit_params (below) will be "locked" in the
%    fitting process.
%    """

function [fitparms] = fit2dgauss (data, params)

% def fitgaussian(data, params, fixed=None, maxfev=0):
%     fixed = fixed or {}

    % Collect necessary values from parameter dict; only those which aren't
    % fixed.
    initial = []
    for param in FIT_PARAMS:
        if param not in fixed:
            if hasattr(params[param], "value"):
                initial.append(params[param].value)
            else:
                initial.append(params[param])

	% NOTE: Matlab subfunction: Can be called only from within this function!
    function residuals(paramlist):
        % Error function to be used in chi-squared fitting
        %	argument paramlist: fitting parameters
        %	type paramlist: numpy.ndarray
        %	argument fixed: parameters to be held frozen
        %	type fixed: dict

        %	returns: 2d-array of diffs. between estimated Gaussian function
        %   and the actual data
        paramlist = list(paramlist)
        gaussian_args = []
        for param in FIT_PARAMS:
            if param in fixed:
                gaussian_args.append(fixed[param])
            else:
                gaussian_args.append(paramlist.pop(0))

        % gaussian() returns a function which takes arguments x, y and returns
        % a Gaussian with parameters gaussian_args evaluated at that point.
        g = gaussian(*gaussian_args)

        % The .compressed() below is essential so the Gaussian fit will not
        % take account of the masked values (=below threshold) at the edges
        % and corners of data (=(masked) array, so rectangular in shape).
		return (numpy.fromfunction(g, data.shape) - data).compressed()

    % maxfev=0, the default, corresponds to 200*(N+1) (NB, not 100*(N+1) as
    % the scipy docs state!) function evaluations, where N is the number of
    % parameters in the solution.
    % Convergence tolerances xtol and ftol established by experiment on images
    % from Paul Hancock's simulations.
    soln, success = scipy.optimize.leastsq(
        residuals, initial, maxfev=maxfev, xtol=1e-4, ftol=1e-4
    )

    if success > 4:
        raise ValueError("leastsq returned %d; bailing out" % (success,))

    % soln contains only the variable parameters; we need to merge the
    % contents of fixed into the soln list.
    % leastsq() returns either a numpy.float64 (if fitting a single value) or
    % a numpy.ndarray (if fitting multiple values); we need to turn that into
    % a list for the merger.
    try:
        % If an ndarray (or other iterable)
        soln = list(soln)
    except TypeError:
        soln = [soln]
    results = fixed.copy()
    for param in FIT_PARAMS:
        if param not in results:
            results[param] = soln.pop(0)

    if results['semiminor'] > results['semimajor']:
        % Swapped axis order is a perfectly valid fit, but inconvenient for
        % the rest of our codebase.
        results['semimajor'], results['semiminor'] = results['semiminor'], results['semimajor']
        results['theta'] += numpy.pi/2

    % Negative axes are a valid fit, since they are squared in the definition
    % of the Gaussian.
    results['semimajor'] = abs(results['semimajor'])
    results['semiminor'] = abs(results['semiminor'])

    return results
