% Function interprets the binary header and returns meta information.
% Arguments:
%	hdr : array of 512 uint8, corresponding to the first 512 bytes of
%		  the visibility record.
% Returns :
%	meta: A structure containing metadata as extracted from the header.
function meta = interpretHdr (obj, hdr)
	% Set the object magic value
    tmp = typecast (hdr(1), 'uint32'); 
	meta.magic  = tmp(1);
    
    if (meta.magic == 999878658)
        meta.tstart = unixtime2mjdsec (hdr(2));
        meta.tend   = unixtime2mjdsec (hdr(3));
    elseif (meta.magic == 0x41415254)
        meta.tstart = unixtime2mjdsec (hdr(3));
        meta.tend   = unixtime2mjdsec (hdr(4));
    end;
end;
