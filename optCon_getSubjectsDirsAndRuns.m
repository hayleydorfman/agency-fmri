function [ subjdirs, nRuns, goodRuns, goodSubjects, subj_original_indices] = optCon_getSubjectsDirsAndRuns()

% use this script for all subj info

% the names of the CORRESPONDING directories from CBS central
subjdirs = {'180508_UCR002','180508_UCR003','180508_UCR004','180508_UCR005','180514_UCR006','180514_UCR007', '180517_UCR008','180517_UCR009','180517_UCR010','180517_UCR011', '180521_UCR012','180521_UCR013', '180522_UCR014', '180523_UCR015', '180523_UCR016', '180523_UCR017', '180524_UCR018', '180530_UCR020','180530_UCR021','180531_UCR022', '180531_UCR023', '180601_UCR024','180601_UCR025', '180604_UCR026', '180605_UCR027','180606_UCR028', '180606_UCR029','180806_UCR030','180813_UCR031', '180814_UCR032', '180815_UCR033', '180816_UCR034', '180817_UCR035', '180820_UCR036', '180822_UCR037', '180822_UCR038', '180823_UCR039', '180828_UCR040'};

subj_original_indices = [2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40];

% assumes runs are always in order: 1,2,3,4,...
nRuns = {}; % runs per subject [We don't need this]


%goodSubjects = [1:4, 5, 6, 8, 9, 10, 11, 13, 15, 16, 17, 18:21, 23:26, 28:30, 32:37]; %this is SPM index!
% this version is excluding the following subjs (in original indices): 8, 13,15,24,29,33,40
% this version was used for initial imaging analyses. It excludes subjs for
% missing trials, missing runs, bad img data quality, and modeling issues
% (i.e., psi parameter always 1 because they only chose one option for
% every trial)
% S1

%goodSubjects = [1, 2, 4, 5, 6, 8, 9, 10, 11, 13, 15, 16, 18, 20, 21, 23, 26, 28:30, 32:36]; %this is SPM index!
% this version is excluding the following subjs (in original indices): 4, 8, 13,15,18,21,24,25,26,29,33,39,40
% this version was used for imaging analyses with subjs who also met
% behavioral accuracy criteria - MOST GLMS for final paper will be this
% subset
%n=25
% S2

goodSubjects = [1, 2, 4, 5, 6, 8, 9, 10, 11, 13, 15, 16, 18, 20, 21, 23, 26, 28:30, 32:34]; %this is SPM index!
% this version is excluding the following subjs (in original indices): 4,
% 8, 13,15,18,21,24,25,26,29,33,39,40, 37, 38
% this version was used for imaging analyses with subjs who also met
% behavioral accuracy criteria - USE FOR GLM 16 (or any GLM where you look
% at agency beliefs) - 2 subjs in the subset above NEVER said that the
% latent agent caused wins
% n = 23
%S3

%goodSubjects = [1, 2, 4, 5, 6, 8, 9, 10, 11, 13, 15, 16, 18, 20, 21, 23, 26, 28, 30, 32:35]; %this is SPM index!
% this version was used for GLM 17; assertions failed for 2 additional
% subjs, so this is a new subset of 23 subjs
%S4

%goodSubjects = [1, 2, 4, 5, 6, 8, 9, 10, 11, 13, 15, 16, 18, 21, 23, 26, 28:30, 32:36]; %this is SPM index!
% this version was used for GLM 17; assertions failed for 2 additional
% subjs, so this is a new subset of 23 subjs
%S5




% which runs to include/exclude for each subject (for everyone! - don't
% change this)
goodRuns = {logical([1 1 1 1]), logical([1 1 1 1]), logical([1 1 1 1]), logical([1 1 1 1]), logical([0 1 1 1]), ...
                logical([1 1 1 1]), logical([0 1 1 1]), logical([1 1 1 1]), logical([0 1 1 1]), logical([1 1 1 1]), ...
                logical([1 1 1 1]), logical([1 0 1 1]), logical([1 1 1 1]), logical([1 1 1 1]), logical([1 1 1 1 ]), ...
                logical([1 1 1 1]), logical([1 1 0 1]), logical([1 1 1 1]), logical([1 1 1 1]), logical([1 1 1 1]), ...
                logical([1 1 1 1]), logical([1 1 1 1]), logical([1 1 1 1]), logical([1 1 1 1]), logical([1 1 1 1]), ...
                logical([1 1 1 1]), logical([1 1 1 0]), logical([1 1 1 1]), logical([1 1 1 1]), logical([1 1 1 1]), ...
                logical([1 1 1 1]), logical([1 1 1 1]), logical([1 1 1 1]), logical([1 1 1 1]), logical([1 1 1 1]), ...
                logical([1 1 1 0]), logical([1 1 1 0]), logical([1 1 1 0])};






% 
% function [ subjects, subjdirs, goodRuns, goodSubjects ] = exploration_getSubjectsDirsAndRuns()
% 
% % Get the list of subjects, subject directories, and number of runs for the
% % fMRI GLM code
% %
% 
% %[data, metadata] = load_data(fullfile('data', 'fmri.csv'), true, getGoodSubjects());
% 
% % the participant id as entered in psychopy
% subjects = {'uep001', 'uep002', 'uep003', 'uep004', 'uep005', ...
%             'uep006', 'uep007', 'uep008', 'uep009', 'uep010', ... 
%             'uep011', 'uep012', 'uep013', 'uep014', 'uep015', ...
%             'uep016', 'uep017', 'uep018', 'uep019', 'uep020', ...
%             'uep021', 'uep022', 'uep023', 'uep024', 'uep025', ...
%             'uep026', 'uep027', 'uep028', 'uep029', 'uep030', ...
%             'uep031'};
% 
% % should be identical to the list of subjects in the csv file
% % and in the same order
% % this is a basic assumption for getGoodSubjects() to work
% % we are listing them here explicitly as a sanity check for the csv file
% %
% %assert(mean(strcmp(subjects, unique(data.participant)')) == 1);
% 
% % the names of the CORRESPONDING directories from CBS central
% subjdirs = {'180725_UEP_001', '180727_UEP_002', '180727_UEP_003',  '180730_UEP_004', '180730_UEP_005', ...
%             '180801_UEP_006', '180802_UEP_007', '180803_UEP_008',  '180803_UEP_009', '180804_UEP010',  ...
%             '180804_UEP_011', '180804_UEP_012', '180804_UEP_013',  '180804_UEP_014', '180804_UEP_015', ...
%             '180805_UEP_016', '180805_UEP_017', '180805_UEP_018_2','180805_UEP_019', '180805_UEP_020', ...
%             '180805_UEP_021', '180806_UEP_022', '180806_UEP_023',  '180807_UEP_024', '180807_UEP_025', ...
%             '180807_UEP_026', '180808_UEP_027', '180808_UEP_028',  '180808_UEP_029', '180809_UEP_030', ...
%             '180809_UEP_031'};
% 

% 
% 
% % which subjects are good
% goodSubjects = 1:31;
%  
%assert(numel(subjects) == numel(subjdirs));
%assert(numel(subjects) == numel(goodRuns));
% 
