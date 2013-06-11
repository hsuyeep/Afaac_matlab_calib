% Given a sequence of numbers, this function estimates the mean of the set by 
% iteratively eliminating values above a parameterized threshold.
% pep/04Apr13

function [m, v, sel] = robustmean (x, sig)
	if (length (x) < 2)
		m = x(1);
		v = 0;
	end;

	sel = ones (1, length (x));
	for ind = 1:10
		m1 = mean (x(sel == 1));
		v1 = std (x(sel == 1));
		sel = (abs (x - m1) / v1) < sig;
	end;	

	m = m1; v = v1;
