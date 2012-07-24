function JD = JulianDay(time)

% JD = JulianDay(time)
%
% Computes the Julian day (including decimal fraction) corresponding to the
% requested time.
%
% Parameter
% time : UTC time as date number as returned by DATENUM function
%
% Return value
% JD  : Julian Day number including decimal fraction
%
% This function was written by Sebastiaan van der Tol in 2002 to support
% his ThEA experiments. The conversion formulas were taken from [1].
%
% Rererences
% [1] Peter Duffet-Smith, "Practical Astronomy with Your Calculator",
%     Cambridge University Press, 1979
%
% SJW, 2006

JD = time - datenum(1998, 2, 13, 12,0,0) + 2450858;  
