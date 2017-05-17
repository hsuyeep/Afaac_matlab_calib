% Function to read in multiple timeslices of data, sigmaclip in time
% axis and write the averaged visibility back to the object. For implementing
% time domain averaging in visibilities before imaging.
% Calls readRec() underneath.
% NOTE: The VisRec object no longer corresponds to the raw visibility in
% the associated file after calling this method! Only the trecstart,
% freq, nchan, acm_xx and acm_yy are updated. recbytesize, datfloatsize
% are not updated to allow reading in the next set of raw records.
% Arguments:
%	pol  : The pols to  be read in. Bool array, [XX, XY, YX, YY]	
%	chans: The channel subset to  be read in. Range [1:obj.nchans]	
%   tavg : The number of seconds (modulo integration time) to average over.
% Returns:
%   dat  : Time averaged visibilities, flagged both along the channel and time
%          axis.
function dat = readTimeAvgRec (obj, pol, chans, tavg)
    
    % Read in the next record to determine its timestamp.
    obj.readRec (pol, chans);

    % Determine the end timestamp
    tend = obj.trecstart + tavg;

    % Check!
    assert (tend > obj.trecstart);

    % Internally store records
    xxtslices(1, :) = obj.xx;
    yytslices(1, :) = obj.yy;
    trecavg = obj.trecstart - obj.dt + tavg/2;

    if (obj.deb > 2)
        fprintf (2, 'readTimeAvgRec: output time %.2f', dat.trecstart);
    end;

    % Ingest records till avg. is reached.
    ind = 2;
    while (obj.trecstart + obj.dt <= tend)
        obj.readRec (pol, chans);
        xxtslices (ind, :) = obj.xx;
        yytslices (ind, :) = obj.yy;
        ind = ind + 1;
        if (obj.deb > 2)
            fprintf (2, '%d:%.2f ', ind, obj.trecstart);
        end;
    end;

    if (obj.timeflag)
        % Selection mask
        selxx = ones (size (xxtslices));
        selyy = ones (size (yytslices));
        % ### UNIMPLEMENTED, DONT USE!
        while 1
            m = nanmean (xxtslices);
            s = nanstd  (xxtslices);
            xxtslices(abs(xxtslices-repmat (m, [obj.nbline, 1])) > repmat (obj.flagsig*abs(s), [obj.nbline, 1])) = nan;

            m = nanmean (yytslices);
            s = nanstd  (yytslices);
            yytslices(abs(yytslices-repmat (m, [obj.nbline, 1])) > repmat (obj.flagsig*abs(s), [obj.nbline, 1])) = nan;
        end;
    else
        obj.xx = mean (xxtslices, 1);
        obj.yy = mean (yytslices, 1);
        obj.trecstart = trecavg;
        % TODO: Check what else needs to be updated in order to keep the
        % the parent VisRec object's metadata consistent.
    end;
end;
