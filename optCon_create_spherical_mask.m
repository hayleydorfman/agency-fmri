% wrapper for ccnl_create_spherical_mask.m

% These lines are for sanity checking the conversion - change numbers to
% cor
[V, Y] = ccnl_extract_clusters(optCon_expt(), 6, 'RPE', 0.001, '+', 0.05, 20, 3);
%Y(35,60,31)

r = 10/1.5; % radius (divide by voxel size)


[mask, Vmask, Ymask] = load_mask('masks/mask.nii'); % load the mask

filename = 'masks/left_VS.nii' % change name - this will be output this way

cor = mni2cor([-14,4,-12], Vmask.mat) %change coordinates depending on which region you want the spherical mask to cover

[sphere_mask, sphere_coords, sphere_vol] = ccnl_create_spherical_mask(cor(1), cor(2), cor(3), r, filename, mask);

view_mask(filename)

Y(cor(1),cor(2),cor(3)) % this should be the same as the one in the table of bspmview


