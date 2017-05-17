% Constructor
% Opens the binary file, initializes internal variables.
% If metainfo. cannot be determined from the data stream,
% uses the externally supplied values.
% Arguments:
%   fname   : Filename of visibility binary file
%   info    : Structure containing observational information.
function  obj = VisRec (fname, info)
	assert (exist (fname) == 2);
	obj.fname = fname;
	try
	    obj.fid  = fopen (obj.fname, 'rb');	
	catch ME
     	fprintf (2, '### Error in opening file %s.\n', obj.fname);
	end;

	if (isfield (info, 'nchan' ) == 0) obj.nchan =  63; else obj.nchan = info.nchan; end;
	if (isfield (info, 'npol'  ) == 0) obj.npol  =   2; else obj.npol = info.npol; end;
	if (isfield (info, 'nelem' ) == 0) obj.nelem = 288; else obj.nelem = info.nelem; end;
    obj.nbline = obj.nelem * (obj.nelem+1)/2;
	% if (isfield (info, 'nbline') == 0) obj.nbline= 41616; else obj.nbline = info.nbline; end;
	if (isfield (info, 'freqflag') == 0) obj.freqflag= 0; else obj.freqflag = info.freqflag; end;
	if (isfield (info, 'timeflag') == 0) obj.timeflag= 0; else obj.timeflag = info.timeflag; end;
	if (isfield (info, 'freq') == 0) obj.freq = 0; else obj.freq = info.freq; end;
	if (isfield (info, 'skip'  ) == 0) obj.skip  =   0; else obj.skip = info.skip; end;
	if (isfield (info, 'deb'   ) == 0) obj.deb   =   0; else obj.deb = info.deb; end;
				
    obj.recbytesize = 0;
    obj.datfloatsize = 0;
    obj.tfilestart = 0;
    obj.tfileend = 0;
    obj.trecstart = 0;
    obj.dt = 0;
    obj = class (obj, "VisRec");
    
	% obj.recdat = RecDat(); % Self contained data structure to hold a single records' data.
	% obj.getVisMeta (info); % Initialize internal data based on input data.
end;
