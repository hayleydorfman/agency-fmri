% Pull beta series from a seed region for PPI analyses

% Use: B = getBetaSeries(EXPT,GLM number,subj,event name, mask filename);
% Example: getBetaSeries(optCon_expt(),7,21,'feedback_onset', 'left_VS.nii');

function B = getBetaSeries(EXPT,GLM,subj,run, event_prefix, mask)

prefix = [event_prefix, '_run_', num2str(run)];

% load mask
    [mask_format, mask, Vmask] = get_mask_format_helper(mask);
    assert(strcmp(mask_format, 'mask'), 'Improper mask');


  % load betas 
            modeldir = fullfile(EXPT.modeldir,['model',num2str(GLM)],['subj',num2str(subj)]);
            load(fullfile(modeldir,'SPM.mat'));
            assert(isequal(SPM.Vbeta(1).dim, Vmask.dim), 'Different dimensions between mask and betas');

            which = contains(SPM.xX.name, prefix); % betas for given event
            %which(which) = rsa.which_betas; % of those, only betas for given trials
            cdir = pwd;
            cd(modeldir); % b/c SPM.Vbeta are relative to modeldir
            B = spm_data_read(SPM.Vbeta(which), find(mask));
            cd(cdir);
            
B = mean(B,2);

assert(size(B,1) > 0, 'no betas - likely wrong prefix');
