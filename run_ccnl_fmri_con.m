
 [allSubjects, subjdirs, goodRuns, goodSubjs, subj_original_indices] = optCon_getSubjectsDirsAndRuns();


% ccnl_fmri_con(optCon_expt(), 1, ...
%    {'benevolent', 'adversarial', 'benevolent - adversarial'}, ...
%     goodSubjs);
% 
% 
% ccnl_fmri_con(optCon_expt(), 2, ...
%    {'wins_benevolent - losses_benevolent', 'wins_adversarial - losses_adversarial', ...
%     'wins_benevolent + wins_adversarial - losses_benevolent - losses_adversarial', ... %all wins minus all losses
%     'wins_benevolent - losses_benevolent - wins_adversarial + losses_adversarial', ...
%     'wins_adversarial - losses_adversarial - wins_benevolent + losses_benevolent', ...
%     'wins_adversarial - wins_benevolent', ...
%     'losses_adversarial - losses_benevolent'}, ... %add to this list the
%     wins-losses
%     goodSubjs);
%         
% % 
% ccnl_fmri_con(optCon_expt(), 3, ...
%    {'benevolent', 'adversarial', 'benevolent - adversarial'}, ...
%     goodSubjs);
% 
% ccnl_fmri_con(optCon_expt(), 4, ...
%    {'RPEpsi'}, ...
%     goodSubjs);
% 
% ccnl_fmri_con(optCon_expt(), 5, ...
%    {'RPE', 'psi', 'psi - RPE', 'RPE - psi'}, ...
%     goodSubjs);
% 
% ccnl_fmri_con(optCon_expt(), 6, ...
%    {'RPE'}, ...
%     goodSubjs);
% 
% ccnl_fmri_con(optCon_expt(), 8, ...
%    {'VS', 'psiVS'}, ...
%     goodSubjs);
% 
% ccnl_fmri_con(optCon_expt(), 9, ...
%    {'VS', 'psiVS'}, ...
%     goodSubjs);
% 
% 
% 
% ccnl_fmri_con(optCon_expt(), 11, ...
%    {'psi'}, ...
%     goodSubjs);
% % 
% ccnl_fmri_con(optCon_expt(), 13, ...
%    {'RPE', 'RPEpsi', 'RPEpsi - RPE', 'RPE - RPEpsi'}, ...
%     goodSubjs);
% 
% ccnl_fmri_con(optCon_expt(), 14, ...
%    {'RPE', 'RPEpsi'}, ...
%     goodSubjs);
% 
% 
% ccnl_fmri_con(optCon_expt(), 15, ...
%    {'RPE', 'psi', 'psi - RPE', 'RPE - psi'}, ...
%     goodSubjs);

%need to use S3
% ccnl_fmri_con(optCon_expt(), 16, ...
%    {'wins_yes', 'wins_no', 'losses_yes', 'losses_no', ...
%    'wins_yes - wins_no', 'losses_yes - losses_no', ...
%    'wins_yes + losses_yes - wins_no - losses_no', ... %all yes guess - all no guess
%    'wins_yes + losses_no - wins_no - losses_yes', ... %you cause good things - i cause good things
%    'wins_yes + wins_no - losses_yes - losses_no', ... % all wins - all losses
%    'wins_yes - losses_yes', 'wins_no - losses_no', ...
%    'wins_yes - losses_no', 'wins_no - losses_yes', ...
%    'losses_no - losses_yes', ... %my fault losses - your fault losses
%    'wins_no - wins_yes - losses_yes + losses_no', ... %interaction1
%    'wins_no - losses_no - wins_yes + losses_yes'}, ... %interaction2 
%     goodSubjs);

% % need to use S4
% ccnl_fmri_con(optCon_expt(), 17, ...
%    {'wins_low', 'wins_high', 'losses_low', 'losses_high', ... %high psi = more self agency, low psi = less agency
%    'wins_high - wins_low', 'losses_high - losses_low', ...
%    'wins_high - losses_high', 'wins_low - losses_low', ...
%    'wins_high + losses_high - wins_low - losses_low', ... %all high psi - all low psi; internal agency - external agency
%    'wins_high + wins_low - losses_high - losses_low'}, ... %all wins - all losses
%     goodSubjs);

% 
% ccnl_fmri_con(optCon_expt(), 19, ...
%    {'psi_low', 'psi_high', 'psi_low - psi_high'}, ... 
%     goodSubjs);

% 
% ccnl_fmri_con(optCon_expt(), 20, ...
%    {'RPE_1_minus_psi'}, ... 
%     goodSubjs);
% 
% ccnl_fmri_con(optCon_expt(), 21, ...
%    {'RPE_1_minus_psi', 'RPEpsi', 'RPE_1_minus_psi - RPEpsi'}, ... 
%     goodSubjs);

% ccnl_fmri_con(optCon_expt(), 23, ...
%    {'yes', 'no', 'yes - no', 'no - yes'}, ... 
%     goodSubjs);


% ccnl_fmri_con(optCon_expt(), 24, ...
%    {'MidTemp', 'adversarialxMidTemp - benevolentxMidTemp'}, ... 
%     goodSubjs);

% ccnl_fmri_con(optCon_expt(), 26, ...
%    {'yes', 'no', 'yes - no', 'no - yes'}, ... 
%     goodSubjs);

% ccnl_fmri_con(optCon_expt(), 27, ...
%    {'yes_benevolent - no_benevolent', 'yes_adversarial - no_adversarial', ...
%    'no_benevolent - yes_benevolent', 'no_adversarial - yes_adversarial', ...
%     'yes_benevolent + yes_adversarial - no_benevolent - no_adversarial', ... %all yes minus all no
%     'yes_benevolent - no_benevolent - yes_adversarial + no_adversarial', ... %interaction 1
%     'yes_adversarial - no_adversarial - yes_benevolent + no_benevolent', ... %interaction 2
%     'yes_adversarial - yes_benevolent', ...
%     'no_adversarial - no_benevolent'}, ... 
%     goodSubjs);

% 
% ccnl_fmri_con(optCon_expt(), 29, ...
%    {'all_trials', 'wins', 'adversarial', 'wins_adversarial'}, ... 
%     goodSubjs);

% ccnl_fmri_con(optCon_expt(), 30, ...
%    {'adversarial'}, ... 
%     goodSubjs);
% % 
% ccnl_fmri_con(optCon_expt(), 31, ...
%    {'wins_adversarial'}, ... 
%     goodSubjs);

% ccnl_fmri_con(optCon_expt(), 32, ...
%    {'wins', 'adversarial', 'wins_adversarial'}, ... 
%     goodSubjs);

% ccnl_fmri_con(optCon_expt(), 33, ...
%    {'all_trials', 'wins', 'adversarial'}, ... 
%     goodSubjs);


% ccnl_fmri_con(optCon_expt(), 34, ...
%    {'all_trials', 'wins', 'adversarial', 'wins_adversarial'}, ... 
%     goodSubjs);



% ccnl_fmri_con(optCon_expt(), 35, ...      %%works
%    {'all_trials'}, ... 
%     goodSubjs);

% ccnl_fmri_con(optCon_expt(), 36, ...       %works
%    {'all_trials', 'wins_adversarial'}, ... 
%     goodSubjs);

% ccnl_fmri_con(optCon_expt(), 37, ...       %
%    {'all_trials', 'wins', 'wins_adversarial'}, ... 
%     goodSubjs);

% 
% ccnl_fmri_con(optCon_expt(), 38, ...
%    {'yes', 'no', 'yes - no', 'no - yes'}, ... 
%     goodSubjs);

% ccnl_fmri_con(optCon_expt(), 40, ...       %
%    {'wins_benevolent', 'wins_adversarial', 'wins_adversarial - wins_benevolent', 'wins_benevolent - wins_adversarial'}, ... 
%     goodSubjs);

ccnl_fmri_con(optCon_expt(), 41, ...       %
   {'losses_benevolent', 'losses_adversarial', 'losses_adversarial - losses_benevolent', 'losses_benevolent - losses_adversarial'}, ... 
    goodSubjs);


% ccnl_fmri_con(optCon_expt(), 42, ...       %
%    {'losses_benevolent', 'wins_adversarial', 'wins_adversarial - losses_benevolent', 'losses_benevolent - wins_adversarial'}, ... 
%     goodSubjs);


% ccnl_fmri_con(optCon_expt(), 43, ...       %
%    {'losses_benevolent', 'wins_adversarial', 'wins_adversarial - losses_benevolent', 'losses_benevolent - wins_adversarial'}, ... 
%     goodSubjs);

% need to use S4
% ccnl_fmri_con(optCon_expt(), 17, ...
%    {'wins_low', 'wins_high', 'losses_low', 'losses_high'}, ... %all wins - all losses
%     goodSubjs);

%need to use S2
% ccnl_fmri_con(optCon_expt(), 18, ...
%    {'psi'}, ...
%     goodSubjs);





% ccnl_fmri_con(optCon_expt(),10, ...
%    {'wins_benevolentxpsi - losses_benevolentxpsi', 'wins_adversarialxpsi - losses_adversarialxpsi', ...
%     'wins_benevolentxpsi + wins_adversarialxpsi - losses_benevolentxpsi - losses_adversarialxpsi', ... %all wins minus all losses
%     'wins_benevolentxpsi - losses_benevolentxpsi - wins_adversarialxpsi + losses_adversarialxpsi', ...
%     'wins_adversarialxpsi - losses_adversarialxpsi - wins_benevolentxpsi + losses_benevolentxpsi', ...
%     'wins_adversarialxpsi - wins_benevolentxpsi', ...
%     'losses_adversarialxpsi - losses_benevolentxpsi', 'psi'}, ... 
%     goodSubjs);





% ccnl_fmri_con(optCon_expt(),12, ...
%    {'winsxpsi - lossesxpsi', 'psi'}, ... 
%     goodSubjs);

