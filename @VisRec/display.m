% Script to display the contents of a VisRec octave class.
% pep/16May17
function display (obj)
    first = true;
    fprintf (1, '%s =', inputname (1));
    fprintf (1, '<-- Operating on file %s\n', obj.fname);
    fprintf (1, '  Start time: %s, End Time: %s\n', datestr (mjdsec2datenum(obj.tfilestart)), datestr (mjdsec2datenum(obj.tfileend)));   
    fprintf (1, '  Operating frequency: %f Hz\n', obj.freq);
    fprintf (1, '  Nant: %d, Npol: %d Nchan:  %d Nbline: %d.\n', obj.nelem, obj.npol, obj.nchan, obj.nbline);
    fprintf (1, '  Recbytesize: %d, datfloatsize: %d.\n', obj.recbytesize, obj.datfloatsize);
    fprintf (1, '  Data sizes xx: %d, acm_xx: %d, yy: %d, acm_yy: %d\n', size (obj.xx), size (obj.acm_xx), size (obj.yy), size (obj.acm_yy));
