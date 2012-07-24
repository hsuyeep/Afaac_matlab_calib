% Script to create an empirical noise map from a given image
% Returns a matrix of the same size as image, with each element containing 
% the variance of a region defined by a rectangular region of size 'rad';
% pep/08May12

function [map] = noisemap (skymap, rad)
    if (rad > length (skymap))
        disp ('Region larger than image!'); return;
    end
    map = zeros (length(skymap));
    for xpos = rad+1:length(skymap)-rad;
        for ypos = rad+1:length(skymap)-rad;
            subarr = skymap (xpos-rad:xpos+rad, ypos-rad:ypos+rad);
            map (xpos, ypos) = var(var(subarr));
        end;
    end;
            
