function BW_mask_shrinked = shrink_BW_mask(BW_mask,nlines)

    BW_mask_shrinked = BW_mask;

    for i=1:nlines
        BW_mask_shrinked(1,:) = [];
        BW_mask_shrinked(:,1) = [];
        BW_mask_shrinked(end,:) = [];
        BW_mask_shrinked(:,end) = [];
    end
     
end