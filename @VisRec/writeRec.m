% Function to write out visibilities as records in the format
% generated by the C++ pipeline, to be inter-operable.
%   NOTE : We write out both polarizations, XX followed by YY.
%          readRec() averages over channels, so we don't allow writing
%          out individual channels for now, only integrated vis.
%          This has to be ensured by the caller.
% Arguments:
%   vis  : Visibility array to write to disk. Can be raw or calibrated
%          vis, both are written out as full matrices.
%   pol  : 0 or 1, indicating vis. is for XX or YY polarity.
%   fid  : The id of the binary file to which to write. New
%          records are always appended to the specified file.
%
% Returns: Status, true for a successful write.
%
% The output header structure matches the one created by the c++
% pipeline, and replicated below. From
% aartfaac-calibration-pipeline/src/server/packet.h
%
% struct output_header_t
% {
%  uint64_t magic;                   ///< magic to determine header
%  double start_time;                ///< start time (unix)
%  double end_time;                  ///< end time (unix)
%  int32_t subband;                  ///< lofar subband
%  int32_t num_dipoles;              ///< number of dipoles (288 or 576)
%  int32_t polarization;             ///< XX=0, YY=1
%  int32_t num_channels;             ///< number of channels (<= 64)
%  float ateam_flux[5];              ///< Ateam fluxes (CasA, CygA, Tau,
%  Vir, Sun)
%  std::bitset<5> ateam;             ///< Ateam active
%  std::bitset<64> flagged_channels; ///< bitset of flagged channels (8
%  byte)
%  std::bitset<576> flagged_dipoles; ///< bitset of flagged dipoles (72
%  byte)
%  uint32_t weights[78];             ///< stationweights n*(n+1)/2, n in
%  {6, 12}
%  uint8_t pad[48];                  ///< 512 byte block
% };
function stat = writeRec (obj,vis, pol, fid)

    % Create metadata hdr in correct format
    % hdr.magic       = uint64 (0x4141525446414143);
    hdr.magic       = uint64 (0x4141525400000000, 'b');
    hdr.start_time  = obj.trecstart;
    hdr.end_time    = obj.trecend;
    hdr.subband     = int32 (obj.freq/195312.5);
    hdr.num_dipoles = int32 (obj.nelem);
    hdr.polarization= int32 (pol); % Write out XX first
    hdr.num_channels= int32 (1);
    hdr.ateam_flux  = single ([0,0,0,0,0]); % Not writing out now.
    hdr.rest        = uint8 (zeros (1, 512-15)); % Set all other fields
                                                 % to 0.

    status = 0;
    try
        fwrite (fid, hdr.magic, '2*uint32');
        fwrite (fid, hdr.start_time, 'double');
        fwrite (fid, hdr.end_time, 'double');
        fwrite (fid, hdr.subband, 'uint32');
        fwrite (fid, hdr.num_dipoles, 'uint32');
        fwrite (fid, hdr.polarization, 'uint32');
        fwrite (fid, hdr.num_channels, 'uint32');
        fwrite (fid, hdr.ateam_flux, '5*single');
        fwrite (fid, hdr.rest, '452*uint8');

        fmt = sprintf ('%d*single', length(vis)*2); % for complex
        fwrite (fid, vis, fmt); % Write out float complex.
    catch ME
        fprintf (2, '### VisRec.m:writeRec(): Error in writing to file id %d.\n', fid);
        status = 1;
    end;

    % return the status
    % return status;
end;
