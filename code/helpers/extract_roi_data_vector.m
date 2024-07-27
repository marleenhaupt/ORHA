function roi_data = extract_roi_data_vector(roi, contrast)

    Y = spm_read_vols(spm_vol(roi),1);
    indx = find(Y>0);
    [x,y,z] = ind2sub(size(Y),indx);

    XYZ = [x y z]';

    roi_data = spm_get_data(contrast, XYZ);

end