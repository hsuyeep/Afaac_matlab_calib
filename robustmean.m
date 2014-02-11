% Given a sequence of numbers, this function estimates the mean of the set by 
% iteratively eliminating values above a parameterized threshold.
% pep/04Apr13
% Updated to accept complex inputs.
% pep/18Jun13
% Arguments:
%	x: Working vector of real or complex numbers.
% sig: Threshold above sigma on which to clip values.
% Returns:
%	m: Robust mean, computed iteratively
%	v: Robust variance, computed over the final remainder set
% sel: Selection vector of size x. 0's indicate vector members which 
%	   have been eliminated.

function [m, v, sel] = robustmean (x, sig)
	if (length (x) < 2)
		m = x(1);
		v = 0;
	end;

	sel = ones (length (x), 1);
	nandat = sel;
	nandat (isnan (x) == 1) = 0; % list of nan values in data

	if (isreal(x) == 1)
		for ind = 1:10
			sel = (sel & nandat); 
			m1 = mean (x(sel == 1));
			v1 = std (x(sel == 1));
			sel = (abs (x - m1) / v1) < sig;
		end;	
		m = m1; v = v1;
	else
		selr = sel; seli = sel;
		for ind = 1:10
			sel = selr & seli & nandat; 
			mr1 = mean (real(x(sel == 1)));
			mi1 = mean (imag(x(sel == 1)));

			vr1 = std (real(x(sel == 1)));
			vi1 = std (imag(x(sel == 1)));
			selr = (abs (real(x) - mr1) / vr1) < sig;
			seli = (abs (imag(x) - mi1) / vi1) < sig;
			% We reject if any of real/imag breaks threshold.
		end;	
		m = mr1 + 1*i*mi1;
		v = vr1 + 1*i*vi1;
	end;
		

