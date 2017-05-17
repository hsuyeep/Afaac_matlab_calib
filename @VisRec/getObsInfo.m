
% Return the internal state of the VisRec object.
function [obs] = getObsInfo (obj)
    if (isfield (obj, 'nchan' ) == 0)    obs.nchan     = obj.nchan; end;    
    if (isfield (obj, 'npol'  ) == 0)    obs.npol      = obj.npol; end;     
    if (isfield (obj, 'nelem' ) == 0)    obs.nelem     = obj.nelem; end;    
    if (isfield (obj, 'nbline') == 0)    obs.nbline    = obj.nbline; end;   
    if (isfield (obj, 'freqflag')== 0)   obs.freqflag  = obj.freqflag; end; 
    if (isfield (obj, 'timeflag')== 0)   obs.timeflag  = obj.timeflag; end; 
    if (isfield (obj, 'freq')   == 0)    obs.freq      = obj.freq; end; 
    if (isfield (obj, 'skip')   == 0)    obs.skip      = obj.skip; end; 
end;
