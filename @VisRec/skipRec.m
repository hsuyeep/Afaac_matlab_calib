% Function to skip raw records in the .vis file.
% Assumes that the observation parameters of data are already
% initialized.
% Arguments:
%   rec : Number of records to skip.
%   unit: string ('rec', 'dt', 'mjdsec') Whether to skip records or seconds, or to
%   an absolute time.
% Returns:
%   Void.
function skipRec (obj, rec, unit)
    % If unit is time (sec), determine number of records to skip based
    % on integration time and current time.
    switch lower(unit)
        case 'sec' % Should be in seconds
            recs2move = int32(rec/obj.dt);
            fprintf (1, '<-- Moving %d recs for %d seconds.\n', recs2move, rec)
            stat = fseek (obj.fid, obj.recbytesize*recs2move, 0);
            if (stat < 0)
                throw (MException ('VisRec.m:skipRec()', ferror (obj.fid)));
            end;
        
        case 'mjdsec'
            assert (obj.tfileend > rec);
            assert (obj.tfilestart < rec);
            recs2move = int32((rec-obj.trecstart)/obj.dt);
            fprintf (1, '<-- Moving %d recs (class %s) for %f seconds.\n', recs2move, class (recs2move), rec - obj.trecstart)
            stat = fseek (obj.fid, obj.recbytesize*recs2move, 0);
            if (stat < 0)
                throw (MException ('VisRec.m:skipRec()', ferror (obj.fid)));
            end;
        
        case 'rec' 
            stat = fseek (obj.fid, obj.recbytesize*rec, 0);
            if (stat < 0)
                throw (MException ('VisRec.m:skipRec()', ferror (obj.fid)));
            end;
    
        otherwise
            fprintf (2, '### Unknown skip option %s.\n', unit);
    end;
end;
