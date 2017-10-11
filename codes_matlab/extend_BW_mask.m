function BW_mask_extended = extend_BW_mask(BW_mask,nlines)

    [height_ini,width_ini] = size(BW_mask);
    BW_mask_extended = [false(nlines,nlines) false(nlines,width_ini) false(nlines,nlines); ...
                        false(height_ini,nlines) BW_mask false(height_ini,nlines); ...
                        false(nlines,nlines) false(nlines,width_ini) false(nlines,nlines)];
       
end