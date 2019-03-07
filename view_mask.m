function view_mask( mask )
% check a mask against the mean structural, make sure it looks right
%
% INPUT:
% mask = path of .nii file with mask
%
EXPT = optCon_expt();
P = {fullfile('glmOutput', 'mean.nii'), fullfile(mask)};
spm_check_registration(char(P));