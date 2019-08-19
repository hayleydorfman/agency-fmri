function multi = optCon_create_multi(glmodel, subj, run, save_output)

    % Create multi structure, helper function for creating EXPT in
    % optCon_expt.m
    %
    % USAGE: multi = optCon_create_multi(model,subj,run)
    %
    % INPUTS:
    %   glmodel - positive integer indicating general linear model
    %   subj - integer specifying which subject is being analyzed
    %   run - integer specifying the run
    %
    % OUTPUTS:
    %   multi - a structure with the following fields
    %        .names{i}
    %        .onsets{i}
    %        .durations{i}
    %        optional:
    %        .pmod(i).name
    %        .pmod(i).param
    %        .pmod(i).poly
    %
    % Momchil Tomov, July 2018
    % Hayley Dorfman, Nov 2018

    if nargin < 4 || isempty(save_output)
        save_output = false;
    end

    fprintf('glm %d, spm subj %d, spm run %d\n', glmodel, subj, run)
    
    %load data
    data = readtable('data.csv');
    

    [allSubjects, subjdirs, goodRuns, goodSubjs, subj_original_indices] = optCon_getSubjectsDirsAndRuns();
    
  
    % skip bad runs
    runs = find(goodRuns{subj});
    run = runs(run);
    fprintf('actual run %03d.nii \n', run);
    
    % match your subj IDs from the data.csv to the index IDs for SPM
    subj_original_idx = subj_original_indices(subj);
    fprintf('subj_original_idx UCR%03d \n', subj_original_idx);
    

    % distribution = results_V.Distribution;
    % dispersion = results_V.Dispersion;
    % dispersion_estimate = results_V.DispersionEstimate;
    % num_variables = results_V.NumVariables;
    % num_predictors = results_V.NumPredictors;
    % num_observation = results_V.NumObservation;
    % variables = results_V.Variables;
    % loglik    = results_V.LogLikelihood;
    % dfe = results_V.DFE;
    % sse = results_V.SSE;
    % sst = results_V.SST;
    % ssr = results_V.SSR;
    % coefcov = results_V.CoefficientCovariance;
    

    
    %can index by data.condition, etc
    % Timings must be in seconds!!!
    
    %subset your df based on things you care about:
    which_trials = data.run_num == run & data.participant == subj_original_idx;

    %define what will go in your glms here (don't forget ' !!!!)
     
    trial_conds = data.condition(which_trials)'; %conditions in your run
    save debug.mat;
    condition = trial_conds{1};
    trial_onsets = data.trial_onset(which_trials)'; %initial trial onset
    reaction_onsets = data.trial_offset(which_trials)'; %reaction onset (when they press the button)
    feedback_onsets = data.fb_onset(which_trials)'; %feedback onset
    latent_onsets = data.interv_onset(which_trials)'; %latent intervention onset
    latent_offsets = data.interv_offset(which_trials)'; %latent intervention offset (when they press the button)
    %ITI onset ??? ask if we should add
    wins = strcmp(data.feedback(which_trials),'1')'; %win trials
    losses = strcmp(data.feedback(which_trials), '0')'; %loss trials
    wins_yes = strcmp(data.feedback(which_trials),'1')' & strcmp(data.latent_guess(which_trials), '1')';
    wins_no = strcmp(data.feedback(which_trials),'1')' & strcmp(data.latent_guess(which_trials), '0')';
    losses_yes = strcmp(data.feedback(which_trials),'0')' & strcmp(data.latent_guess(which_trials), '1')';
    losses_no = strcmp(data.feedback(which_trials),'0')' & strcmp(data.latent_guess(which_trials), '0')';
    yes = strcmp(data.latent_guess(which_trials), '1')';
    no = strcmp(data.latent_guess(which_trials), '0')';
    adversarial = strcmp(data.condition(which_trials), 'adversarial')';
    benevolent = strcmp(data.condition(which_trials), 'benevolent')';
    wins_adversarial = strcmp(data.feedback(which_trials),'1')' & strcmp(data.condition(which_trials), 'adversarial')';
    timeouts = strcmp(data.feedback(which_trials), 'NA')'; %missing trials
    
    % timeouts_latent_guess is used because when there is a timeout from
    % feedback, you automatically miss a latent guess, BUT if you dont time
    % out for feedback, you can still ALSO timeout a latent guess
    timeouts_latent_guess = timeouts | strcmp(data.latent_guess(which_trials), 'NA')'; %missing latent guesses
    assert(sum(wins) + sum(losses) + sum(timeouts) == sum(which_trials));

    % GLMs
    %
    switch glmodel

        % Benevolent - Adversarial
        % Condition @ trial_onset 
        % nuisance @ choice and feedback onset 
        % 
        case 1 
            % condition @ trial onset
            %
            
    
            multi.names{1} = condition;
            multi.onsets{1} = trial_onsets;
            multi.durations{1} = zeros(size(multi.onsets{1})); %impulse regressor (as opposed to boxcar regressor)


            % nuisance @ reaction onset
            %
            multi.names{2} = 'reaction_onset';
            multi.onsets{2} = reaction_onsets;
            multi.durations{2} = zeros(size(multi.onsets{2}));

            % nuisance @ feedback onset
            %
            multi.names{3} = 'feedback_onset';
            multi.onsets{3} = feedback_onsets;
            multi.durations{3} = zeros(size(multi.onsets{3}));
            
                        % nuisance @ latent onset
            %
            multi.names{4} = 'latent_onset';
            multi.onsets{4} = latent_onsets;
            multi.durations{4} = zeros(size(multi.onsets{4}));
            
                        % nuisance @ latent offset
            %
            multi.names{5} = 'latent_offset';
            multi.onsets{5} = latent_offsets;
            multi.durations{5} = zeros(size(multi.onsets{5}));
         
         
            
        % Win vs. Loss in Benevolent vs. Adversarial
        % Condition @ feedback onset 
        % nuisance @ reaction and trial onset 
        % Contrast to run: 1) wins_benevolent - losses_benevolent, 2)
        % wins_adversarial - losses_adversarial, 3)wins_benevolent +
        % wins_adversarial - losses_benevolent - losses_adversarial (so all
        % wins - all losses), 4)wins_adversarial - losses_adversarial -
        % wins_benevolent + losses_benevolent
        
        case 2 
            
            event_index = 1;
            
             if any(wins)         
                % wins for condition
                multi.names{event_index} = ['wins_',condition];
                multi.onsets{event_index} = feedback_onsets(wins);
                multi.durations{event_index} = zeros(size(multi.onsets{event_index})); %impulse regressor (as opposed to boxcar regressor)
                event_index = event_index +1;
             end
           
            if any(losses)
                % losses for condition
                multi.names{event_index} = ['losses_',condition];
                multi.onsets{event_index} = feedback_onsets(losses);
                multi.durations{event_index} = zeros(size(multi.onsets{event_index})); %impulse regressor (as opposed to boxcar regressor)
                event_index = event_index +1;
            end

            % nuisance @ reaction onset
            %
            multi.names{event_index} = 'reaction_onset';
            multi.onsets{event_index} = reaction_onsets;
            multi.durations{event_index} = zeros(size(multi.onsets{event_index}));
            event_index = event_index +1;

            % nuisance @ trial onset
            %
            multi.names{event_index} = 'trial_onset';
            multi.onsets{event_index} = trial_onsets;
            multi.durations{event_index} = zeros(size(multi.onsets{event_index}));
            event_index = event_index +1;
            
            % nuisance @ latent onset
            %
            multi.names{event_index} = 'latent_onset';
            multi.onsets{event_index} = latent_onsets;
            multi.durations{event_index} = zeros(size(multi.onsets{event_index}));
            
                        % nuisance @ latent offset
            %
            multi.names{event_index} = 'latent_offset';
            multi.onsets{event_index} = latent_offsets;
            multi.durations{event_index} = zeros(size(multi.onsets{event_index}));

         
               
            if sum(timeouts) > 0
               multi.names{event_index} = 'feedback_onset_timeouts';
               multi.onsets{event_index} = feedback_onsets(timeouts); % timeouts only
               multi.durations{event_index} = zeros(size(multi.onsets{event_index}));
            end
            
            
        % Same as GLM1 but at feedback onset (Benevolent - Adversarial)
        % Condition @ fb_onset 
        % nuisance @ trial and choice onset 
        % 
        case 3 
            % condition @ fb onset
            %
            
    
            multi.names{1} = condition;
            multi.onsets{1} = feedback_onsets;
            multi.durations{1} = zeros(size(multi.onsets{1})); %impulse regressor (as opposed to boxcar regressor)


            % nuisance @ reaction onset
            %
            multi.names{2} = 'reaction_onset';
            multi.onsets{2} = reaction_onsets;
            multi.durations{2} = zeros(size(multi.onsets{2}));

            % nuisance @ trial onset
            %
            multi.names{3} = 'trial_onset';
            multi.onsets{3} = trial_onsets;
            multi.durations{3} = zeros(size(multi.onsets{3}));
            
            % nuisance @ latent onset
            %
            multi.names{4} = 'latent_onset';
            multi.onsets{4} = latent_onsets;
            multi.durations{4} = zeros(size(multi.onsets{4}));
            
            % nuisance @ latent offset
            %
            multi.names{5} = 'latent_offset';
            multi.onsets{5} = latent_offsets;
            multi.durations{5} = zeros(size(multi.onsets{5}));

           
            
        % RPE * PSI (multiplicative) @ feedback onset
        % nuisance @ trial and choice onset 
        % 
        case 4 
            % 
            %
            
    
            multi.names{1} = condition;
            multi.onsets{1} = feedback_onsets(~timeouts_latent_guess);
            multi.durations{1} = zeros(size(multi.onsets{1})); %impulse regressor (as opposed to boxcar regressor)

           multi.orth{1} = 0; % do not orthogonalise them  
           
           [RPE,psi] = get_latents(subj_original_idx, run);
           
           multi.pmod(1).name{1} = 'RPEpsi';
           multi.pmod(1).param{1} = RPE' .* psi';
           multi.pmod(1).poly{1} = 1;   
           

            % nuisance @ reaction onset
            %
            multi.names{2} = 'reaction_onset';
            multi.onsets{2} = reaction_onsets;
            multi.durations{2} = zeros(size(multi.onsets{2}));

            % nuisance @ trial onset
            %
            multi.names{3} = 'trial_onset';
            multi.onsets{3} = trial_onsets;
            multi.durations{3} = zeros(size(multi.onsets{3}));
            
            % nuisance @ latent onset
            %
            multi.names{4} = 'latent_onset';
            multi.onsets{4} = latent_onsets;
            multi.durations{4} = zeros(size(multi.onsets{4}));
            
            % nuisance @ latent offset
            %
            multi.names{5} = 'latent_offset';
            multi.onsets{5} = latent_offsets;
            multi.durations{5} = zeros(size(multi.onsets{5}));

            
               
            if sum(timeouts_latent_guess) > 0
               multi.names{6} = 'feedback_onset_timeouts';
               multi.onsets{6} = feedback_onsets(timeouts_latent_guess); % timeouts only
               multi.durations{6} = zeros(size(multi.onsets{6}));
            end     
            
            
            
            
            
       % RPE and PSI (additive) @ feedback onset
        % nuisance @ trial and choice onset 
        % 
        case 5 
            % 
            %
            
    
            multi.names{1} = condition;
            multi.onsets{1} = feedback_onsets(~timeouts_latent_guess);
            multi.durations{1} = zeros(size(multi.onsets{1})); %impulse regressor (as opposed to boxcar regressor)

           multi.orth{1} = 0; % do not orthogonalise them  
           
           [RPE,psi] = get_latents(subj_original_idx, run);
                  
          
           multi.pmod(1).name{1} = 'psi';
           multi.pmod(1).param{1} = psi';
           multi.pmod(1).poly{1} = 1; 

           multi.pmod(1).name{2} = 'RPE';
           multi.pmod(1).param{2} = RPE';
           multi.pmod(1).poly{2} = 1; 

            % nuisance @ reaction onset
            %
            multi.names{2} = 'reaction_onset';
            multi.onsets{2} = reaction_onsets;
            multi.durations{2} = zeros(size(multi.onsets{2}));

            % nuisance @ feedback onset
            %
            multi.names{3} = 'trial_onset';
            multi.onsets{3} = trial_onsets;
            multi.durations{3} = zeros(size(multi.onsets{3}));
            
                        % nuisance @ latent onset
            %
            multi.names{4} = 'latent_onset';
            multi.onsets{4} = latent_onsets;
            multi.durations{4} = zeros(size(multi.onsets{4}));
            
            % nuisance @ latent offset
            %
            multi.names{5} = 'latent_offset';
            multi.onsets{5} = latent_offsets;
            multi.durations{5} = zeros(size(multi.onsets{5}));
            
               
            if sum(timeouts_latent_guess) > 0
               multi.names{6} = 'feedback_onset_timeouts';
               multi.onsets{6} = feedback_onsets(timeouts_latent_guess); % timeouts only
               multi.durations{6} = zeros(size(multi.onsets{6}));
            end  
            
            
            
             % RPE only @ feedback onset
        % nuisance @ trial and choice onset 
        % 
        case 6
            % 
            %
            
    
            multi.names{1} = condition;
            multi.onsets{1} = feedback_onsets(~timeouts_latent_guess);
            multi.durations{1} = zeros(size(multi.onsets{1})); %impulse regressor (as opposed to boxcar regressor)

           multi.orth{1} = 0; % do not orthogonalise them  
           
           [RPE,psi] = get_latents(subj_original_idx, run);

           multi.pmod(1).name{1} = 'RPE';
           multi.pmod(1).param{1} = RPE';
           multi.pmod(1).poly{1} = 1; 

            % nuisance @ reaction onset
            %
            multi.names{2} = 'reaction_onset';
            multi.onsets{2} = reaction_onsets;
            multi.durations{2} = zeros(size(multi.onsets{2}));

            % nuisance @ feedback onset
            %
            multi.names{3} = 'trial_onset';
            multi.onsets{3} = trial_onsets;
            multi.durations{3} = zeros(size(multi.onsets{3}));
            
            % nuisance @ latent onset
            %
            multi.names{4} = 'latent_onset';
            multi.onsets{4} = latent_onsets;
            multi.durations{4} = zeros(size(multi.onsets{4}));
            
            % nuisance @ latent offset
            %
            multi.names{5} = 'latent_offset';
            multi.onsets{5} = latent_offsets;
            multi.durations{5} = zeros(size(multi.onsets{5}));
            
               
            if sum(timeouts_latent_guess) > 0
               multi.names{6} = 'feedback_onset_timeouts';
               multi.onsets{6} = feedback_onsets(timeouts_latent_guess); % timeouts only
               multi.durations{6} = zeros(size(multi.onsets{6}));
            end  
            
            
         % Beta series extraction for PPI (part 1/2)
        %  
        % 
        case 7
            % 
            %
           idx = 0;
           trial = data.trial(which_trials);
           
            for t = 1:numel(feedback_onsets)
               idx = idx + 1;
               suffix = ['run_', num2str(run), '_trial_', num2str(trial(t))];
               multi.names{idx} = ['feedback_onset_', suffix];
               multi.onsets{idx} = [feedback_onsets(t)];
               multi.durations{idx} = [0];
            end
           
%            for t = 1:numel(reaction_onsets)
%                idx = idx + 1;
%                suffix = ['run_', num2str(run), '_trial_', num2str(trial(t))];
%                multi.names{idx} = ['reaction_onset_', suffix];
%                multi.onsets{idx} = [reaction_onsets(t)];
%                multi.durations{idx} = [0];
%            end
           
          for t = 1:numel(trial_onsets)
               idx = idx + 1;
               suffix = ['run_', num2str(run), '_trial_', num2str(trial(t))];
               multi.names{idx} = ['trial_onset_', suffix];
               multi.onsets{idx} = [trial_onsets(t)];
               multi.durations{idx} = [0];
          end
      
         
       % Actual PPI (part 2/2) - for left VS 
       % Test if the degree of belief of agency (psi from model) modulates
       % strength of correlation between seed region & every other region
       % in brain
       % seed region = striatum (from previous GLM)
        case 8    
            
            multi.names{1} = condition;
            multi.onsets{1} = feedback_onsets(~timeouts_latent_guess);
            multi.durations{1} = zeros(size(multi.onsets{1})); %impulse regressor (as opposed to boxcar regressor)

           multi.orth{1} = 0; % do not orthogonalise them  
           
           % see pg 134 in FMRI book for info on PPI
           % PPI model is: Y=B0+BX+BR+BR*X+E
           
           [RPE,psi] = get_latents(subj_original_idx, run);
           %disp(psi)
           
           VS = getBetaSeries(optCon_expt(),7,subj, run, 'feedback_onset', 'masks/left_VS.nii');
           VS = VS(~timeouts_latent_guess);
        
           multi.pmod(1).name{1} = 'psi'; % this is the X
           multi.pmod(1).param{1} = psi';
           multi.pmod(1).poly{1} = 1;
           
           multi.pmod(1).name{2} = 'VS'; %this is the R
           multi.pmod(1).param{2} = VS';
           multi.pmod(1).poly{2} = 1;   
           
           multi.pmod(1).name{3} = 'psiVS'; %this is R*X
           %disp(psi);
           %disp(VS);
           %save('debug.mat');
           multi.pmod(1).param{3} = VS' .* psi';
           multi.pmod(1).poly{3} = 1;   
           

            % nuisance @ reaction onset
            %
            multi.names{2} = 'reaction_onset';
            multi.onsets{2} = reaction_onsets;
            multi.durations{2} = zeros(size(multi.onsets{2}));

            % nuisance @ feedback onset
            %
            multi.names{3} = 'trial_onset';
            multi.onsets{3} = trial_onsets;
            multi.durations{3} = zeros(size(multi.onsets{3}));
            
            % nuisance @ latent onset
            %
            multi.names{4} = 'latent_onset';
            multi.onsets{4} = latent_onsets;
            multi.durations{4} = zeros(size(multi.onsets{4}));
            
            % nuisance @ latent offset
            %
            multi.names{5} = 'latent_offset';
            multi.onsets{5} = latent_offsets;
            multi.durations{5} = zeros(size(multi.onsets{5}));
            
               
            if sum(timeouts_latent_guess) > 0
               multi.names{6} = 'feedback_onset_timeouts';
               multi.onsets{6} = feedback_onsets(timeouts_latent_guess); % timeouts only
               multi.durations{6} = zeros(size(multi.onsets{6}));
            end     
            
            
        
            
      % Actual PPI (part 2/2) - for right VS 
       % Test if the degree of belief of agency (psi from model) modulates
       % strength of correlation between seed region & every other region
       % in brain
       % seed region = striatum (from previous GLM)
        case 9    
            
            multi.names{1} = condition;
            multi.onsets{1} = feedback_onsets(~timeouts_latent_guess);
            multi.durations{1} = zeros(size(multi.onsets{1})); %impulse regressor (as opposed to boxcar regressor)

           multi.orth{1} = 0; % do not orthogonalise them  
           
           % see pg 134 in FMRI book for info on PPI
           % PPI model is: Y=B0+BX+BR+BR*X+E
           
           [RPE,psi] = get_latents(subj_original_idx, run);
           VS = getBetaSeries(optCon_expt(),7,subj, run, 'feedback_onset', 'masks/right_VS.nii'); 
           VS = VS(~timeouts_latent_guess);
           
           multi.pmod(1).name{1} = 'psi'; % this is the X
           multi.pmod(1).param{1} = psi';
           multi.pmod(1).poly{1} = 1;
           
           multi.pmod(1).name{2} = 'VS'; %this is the R
           multi.pmod(1).param{2} = VS';
           multi.pmod(1).poly{2} = 1;   
           
           multi.pmod(1).name{3} = 'psiVS'; %this is R*X
           multi.pmod(1).param{3} = VS' .* psi';
           multi.pmod(1).poly{3} = 1;   
           

            % nuisance @ reaction onset
            %
            multi.names{2} = 'reaction_onset';
            multi.onsets{2} = reaction_onsets;
            multi.durations{2} = zeros(size(multi.onsets{2}));

            % nuisance @ feedback onset
            %
            multi.names{3} = 'trial_onset';
            multi.onsets{3} = trial_onsets;
            multi.durations{3} = zeros(size(multi.onsets{3}));
            
            % nuisance @ latent onset
            %
            multi.names{4} = 'latent_onset';
            multi.onsets{4} = latent_onsets;
            multi.durations{4} = zeros(size(multi.onsets{4}));
            
            % nuisance @ latent offset
            %
            multi.names{5} = 'latent_offset';
            multi.onsets{5} = latent_offsets;
            multi.durations{5} = zeros(size(multi.onsets{5}));
            
               
            if sum(timeouts_latent_guess) > 0
               multi.names{6} = 'feedback_onset_timeouts';
               multi.onsets{6} = feedback_onsets(timeouts_latent_guess); % timeouts only
               multi.durations{6} = zeros(size(multi.onsets{6}));
            end         
            
            

         % Psi only @ feedback onset
        % nuisance @ trial and choice onset 
        % 
        case 11
            % 
            %
            
    
            multi.names{1} = condition;
            multi.onsets{1} = feedback_onsets(~timeouts_latent_guess);
            multi.durations{1} = zeros(size(multi.onsets{1})); %impulse regressor (as opposed to boxcar regressor)

           multi.orth{1} = 0; % do not orthogonalise them  
           
           [RPE,psi] = get_latents(subj_original_idx, run);

           multi.pmod(1).name{1} = 'psi';
           multi.pmod(1).param{1} = psi';
           multi.pmod(1).poly{1} = 1; 

            % nuisance @ reaction onset
            %
            multi.names{2} = 'reaction_onset';
            multi.onsets{2} = reaction_onsets;
            multi.durations{2} = zeros(size(multi.onsets{2}));

            % nuisance @ feedback onset
            %
            multi.names{3} = 'trial_onset';
            multi.onsets{3} = trial_onsets;
            multi.durations{3} = zeros(size(multi.onsets{3}));
            
            % nuisance @ latent onset
            %
            multi.names{4} = 'latent_onset';
            multi.onsets{4} = latent_onsets;
            multi.durations{4} = zeros(size(multi.onsets{4}));
            
            % nuisance @ latent offset
            %
            multi.names{5} = 'latent_offset';
            multi.onsets{5} = latent_offsets;
            multi.durations{5} = zeros(size(multi.onsets{5}));
            
               
            if sum(timeouts_latent_guess) > 0
               multi.names{6} = 'feedback_onset_timeouts';
               multi.onsets{6} = feedback_onsets(timeouts_latent_guess); % timeouts only
               multi.durations{6} = zeros(size(multi.onsets{6}));
            end
            
            
            
        % RPE * PSI (multiplicative) and RPE @ feedback onset
        % nuisance @ trial and choice onset 
        % 
        case 13 
            % 
            %
            
    
            multi.names{1} = condition;
            multi.onsets{1} = feedback_onsets(~timeouts_latent_guess);
            multi.durations{1} = zeros(size(multi.onsets{1})); %impulse regressor (as opposed to boxcar regressor)

           multi.orth{1} = 0; % do not orthogonalise them  
           
           [RPE,psi] = get_latents(subj_original_idx, run);
           
           multi.pmod(1).name{1} = 'RPEpsi';
           multi.pmod(1).param{1} = RPE' .* psi';
           multi.pmod(1).poly{1} = 1; 
           
           multi.pmod(1).name{2} = 'RPE';
           multi.pmod(1).param{2} = RPE';
           multi.pmod(1).poly{2} = 1;   
           

            % nuisance @ reaction onset
            %
            multi.names{2} = 'reaction_onset';
            multi.onsets{2} = reaction_onsets;
            multi.durations{2} = zeros(size(multi.onsets{2}));

            % nuisance @ feedback onset
            %
            multi.names{3} = 'trial_onset';
            multi.onsets{3} = trial_onsets;
            multi.durations{3} = zeros(size(multi.onsets{3}));
            
            % nuisance @ latent onset
            %
            multi.names{4} = 'latent_onset';
            multi.onsets{4} = latent_onsets;
            multi.durations{4} = zeros(size(multi.onsets{4}));
            
            % nuisance @ latent offset
            %
            multi.names{5} = 'latent_offset';
            multi.onsets{5} = latent_offsets;
            multi.durations{5} = zeros(size(multi.onsets{5}));
            
               
            if sum(timeouts_latent_guess) > 0
               multi.names{6} = 'feedback_onset_timeouts';
               multi.onsets{6} = feedback_onsets(timeouts_latent_guess); % timeouts only
               multi.durations{6} = zeros(size(multi.onsets{6}));
            end     
            
            
            
        % RPE * PSI (multiplicative) and RPE @ feedback onset
        % nuisance @ trial and choice onset 
        % 
        case 14 
            % 
            %
            
    
            multi.names{1} = condition;
            multi.onsets{1} = feedback_onsets(~timeouts_latent_guess);
            multi.durations{1} = zeros(size(multi.onsets{1})); %impulse regressor (as opposed to boxcar regressor)

           multi.orth{1} = 1; % orthogonalise them  
           
           [RPE,psi] = get_latents(subj_original_idx, run);

           
           multi.pmod(1).name{1} = 'RPE';
           multi.pmod(1).param{1} = RPE';
           multi.pmod(1).poly{1} = 1;   
           
           multi.pmod(1).name{2} = 'RPEpsi';
           multi.pmod(1).param{2} = RPE' .* psi';
           multi.pmod(1).poly{2} = 1;
           

            % nuisance @ reaction onset
            %
            multi.names{2} = 'reaction_onset';
            multi.onsets{2} = reaction_onsets;
            multi.durations{2} = zeros(size(multi.onsets{2}));

            % nuisance @ feedback onset
            %
            multi.names{3} = 'trial_onset';
            multi.onsets{3} = trial_onsets;
            multi.durations{3} = zeros(size(multi.onsets{3}));
            
            % nuisance @ latent onset
            %
            multi.names{4} = 'latent_onset';
            multi.onsets{4} = latent_onsets;
            multi.durations{4} = zeros(size(multi.onsets{4}));
            
            % nuisance @ latent offset
            %
            multi.names{5} = 'latent_offset';
            multi.onsets{5} = latent_offsets;
            multi.durations{5} = zeros(size(multi.onsets{5}));
            
               
            if sum(timeouts_latent_guess) > 0
               multi.names{6} = 'feedback_onset_timeouts';
               multi.onsets{6} = feedback_onsets(timeouts_latent_guess); % timeouts only
               multi.durations{6} = zeros(size(multi.onsets{6}));
            end     
   
            
            
            % RPE and PSI (additive) @ feedback onset, now with
            % orthogonalization
        % nuisance @ trial and choice onset 
        % 
        case 15 
            % 
            %
            
    
            multi.names{1} = condition;
            multi.onsets{1} = feedback_onsets(~timeouts_latent_guess);
            multi.durations{1} = zeros(size(multi.onsets{1})); %impulse regressor (as opposed to boxcar regressor)

           multi.orth{1} = 1; % orthogonalise them  
           
           [RPE,psi] = get_latents(subj_original_idx, run);
           
           
           multi.pmod(1).name{1} = 'RPE';
           multi.pmod(1).param{1} = RPE';
           multi.pmod(1).poly{1} = 1; 
          
           multi.pmod(1).name{2} = 'psi';
           multi.pmod(1).param{2} = psi';
           multi.pmod(1).poly{2} = 1; 



            % nuisance @ reaction onset
            %
            multi.names{2} = 'reaction_onset';
            multi.onsets{2} = reaction_onsets;
            multi.durations{2} = zeros(size(multi.onsets{2}));

            % nuisance @ feedback onset
            %
            multi.names{3} = 'trial_onset';
            multi.onsets{3} = trial_onsets;
            multi.durations{3} = zeros(size(multi.onsets{3}));
            
            % nuisance @ latent onset
            %
            multi.names{4} = 'latent_onset';
            multi.onsets{4} = latent_onsets;
            multi.durations{4} = zeros(size(multi.onsets{4}));
            
            % nuisance @ latent offset
            %
            multi.names{5} = 'latent_offset';
            multi.onsets{5} = latent_offsets;
            multi.durations{5} = zeros(size(multi.onsets{5}));
            
               
            if sum(timeouts_latent_guess) > 0
               multi.names{6} = 'feedback_onset_timeouts';
               multi.onsets{6} = feedback_onsets(timeouts_latent_guess); % timeouts only
               multi.durations{6} = zeros(size(multi.onsets{6}));
            end       

            
                
        % Win vs. Loss in Yes vs No
        % Latent guess @ feedback onset 
        % nuisance @ reaction and trial onset and latent
        % offset(button press)

        
        case 16 
            
            event_index = 1;
            
            save debug.mat; 
            
            assert(sum(wins_yes) + sum(wins_no) + sum(losses_yes) + sum(losses_no) <= length(feedback_onsets))
            assert(length(wins_yes) == length(feedback_onsets))
            assert(length(wins_no) == length(feedback_onsets))
            assert(length(losses_yes) == length(feedback_onsets))
            assert(length(losses_no) == length(feedback_onsets))
            
             if any(wins_yes)         
                % wins for latent guess = yes
                multi.names{event_index} = ['wins_yes'];
                multi.onsets{event_index} = feedback_onsets(wins_yes);
                multi.durations{event_index} = zeros(size(multi.onsets{event_index})); %impulse regressor (as opposed to boxcar regressor)
                event_index = event_index +1;
             end
             
             if any(wins_no)         
                % wins for latent guess = no
                multi.names{event_index} = ['wins_no'];
                multi.onsets{event_index} = feedback_onsets(wins_no);
                multi.durations{event_index} = zeros(size(multi.onsets{event_index})); %impulse regressor (as opposed to boxcar regressor)
                event_index = event_index +1;
             end
             
             if any(losses_yes)         
                % losses for latent guess = yes
                multi.names{event_index} = ['losses_yes'];
                multi.onsets{event_index} = feedback_onsets(losses_yes);
                multi.durations{event_index} = zeros(size(multi.onsets{event_index})); %impulse regressor (as opposed to boxcar regressor)
                event_index = event_index +1;
             end
             
             if any(losses_no)         
                % losses for latent guess = no
                multi.names{event_index} = ['losses_no'];
                multi.onsets{event_index} = feedback_onsets(losses_no);
                multi.durations{event_index} = zeros(size(multi.onsets{event_index})); %impulse regressor (as opposed to boxcar regressor)
                event_index = event_index +1;
                
             end
             

            % nuisance @ reaction onset
            %
            multi.names{event_index} = 'reaction_onset';
            multi.onsets{event_index} = reaction_onsets;
            multi.durations{event_index} = zeros(size(multi.onsets{event_index}));
            event_index = event_index +1;

            % nuisance @ trial onset
            %
            multi.names{event_index} = 'trial_onset';
            multi.onsets{event_index} = trial_onsets;
            multi.durations{event_index} = zeros(size(multi.onsets{event_index}));
            event_index = event_index +1;
            
            % nuisance @ latent onset
            %
            multi.names{event_index} = 'latent_onset';
            multi.onsets{event_index} = latent_onsets;
            multi.durations{event_index} = zeros(size(multi.onsets{event_index}));
            
            % nuisance @ latent offset
            %
            multi.names{event_index} = 'latent_offset';
            multi.onsets{event_index} = latent_offsets;
            multi.durations{event_index} = zeros(size(multi.onsets{event_index}));
         
               
            if sum(timeouts) > 0
               multi.names{event_index} = 'feedback_onset_timeouts';
               multi.onsets{event_index} = feedback_onsets(timeouts); % timeouts only
               multi.durations{event_index} = zeros(size(multi.onsets{event_index}));
            end
            
     
            
         % Win vs. Loss in Yes vs No
        % Condition @ feedback onset 
        % nuisance @ reaction and trial onset 
        % Contrast to run: 1) wins_benevolent - losses_benevolent, 2)
        % wins_adversarial - losses_adversarial, 3)wins_benevolent +
        % wins_adversarial - losses_benevolent - losses_adversarial (so all
        % wins - all losses), 4)wins_adversarial - losses_adversarial -
        % wins_benevolent + losses_benevolent     
       
           case 17 
            
            event_index = 1;
            
            [RPE,psi] = get_latents(subj_original_idx, run);
            
            psi_low = psi < 0.52;
            psi_high = ~psi_low;
            
            psi_low_1 = logical(zeros(size(feedback_onsets)));
            psi_low_1(~timeouts_latent_guess) = psi_low;
            
            psi_high_1 = logical(zeros(size(feedback_onsets)));
            psi_high_1(~timeouts_latent_guess) = psi_high;
            
            wins_low = wins & psi_low_1;
            wins_high = wins & psi_high_1;
            losses_low = losses & psi_low_1;
            losses_high = losses & psi_high_1;
                    
            save debug.mat
            
            assert(sum(wins_low) + sum(wins_high) + sum(losses_low) + sum(losses_high) <= length(feedback_onsets))
            assert(length(wins_low) == length(feedback_onsets))
            assert(length(wins_high) == length(feedback_onsets))
            assert(length(losses_low) == length(feedback_onsets))
            assert(length(losses_high) == length(feedback_onsets))
            
             if any(wins_low)         
                % wins for psi low
                multi.names{event_index} = ['wins_low'];
                multi.onsets{event_index} = feedback_onsets(wins_low);
                multi.durations{event_index} = zeros(size(multi.onsets{event_index})); %impulse regressor (as opposed to boxcar regressor)
                event_index = event_index +1;
             end
             
             if any(wins_high)         
                % wins for psi high
                multi.names{event_index} = ['wins_high'];
                multi.onsets{event_index} = feedback_onsets(wins_high);
                multi.durations{event_index} = zeros(size(multi.onsets{event_index})); %impulse regressor (as opposed to boxcar regressor)
                event_index = event_index +1;
             end
             
            if any(losses_low)         
                % losses for psi low
                multi.names{event_index} = ['losses_low'];
                multi.onsets{event_index} = feedback_onsets(losses_low);
                multi.durations{event_index} = zeros(size(multi.onsets{event_index})); %impulse regressor (as opposed to boxcar regressor)
                event_index = event_index +1;
             end
             
             if any(losses_high)         
                % wins for psi high
                multi.names{event_index} = ['losses_high'];
                multi.onsets{event_index} = feedback_onsets(losses_high);
                multi.durations{event_index} = zeros(size(multi.onsets{event_index})); %impulse regressor (as opposed to boxcar regressor)
                event_index = event_index +1;
             end         
             
          

            % nuisance @ reaction onset
            %
            multi.names{event_index} = 'reaction_onset';
            multi.onsets{event_index} = reaction_onsets;
            multi.durations{event_index} = zeros(size(multi.onsets{event_index}));
            event_index = event_index +1;

            % nuisance @ trial onset
            %
            multi.names{event_index} = 'trial_onset';
            multi.onsets{event_index} = trial_onsets;
            multi.durations{event_index} = zeros(size(multi.onsets{event_index}));
            event_index = event_index +1;

           % nuisance @ latent onset
            %
            multi.names{event_index} = 'latent_onset';
            multi.onsets{event_index} = latent_onsets;
            multi.durations{event_index} = zeros(size(multi.onsets{event_index}));
            
            % nuisance @ latent offset
            %
            multi.names{event_index} = 'latent_offset';
            multi.onsets{event_index} = latent_offsets;
            multi.durations{event_index} = zeros(size(multi.onsets{event_index}));
         
               
            if sum(timeouts) > 0
               multi.names{event_index} = 'feedback_onset_timeouts';
               multi.onsets{event_index} = feedback_onsets(timeouts); % timeouts only
               multi.durations{event_index} = zeros(size(multi.onsets{event_index}));
            end
            
               
        % Psi only @ latent onset
        % nuisance @ trial and choice and feedback onset 
        % 
        case 18
            % 
            %
            
    
            multi.names{1} = condition;
            multi.onsets{1} = latent_onsets(~timeouts_latent_guess);
            multi.durations{1} = zeros(size(multi.onsets{1})); %impulse regressor (as opposed to boxcar regressor)

           multi.orth{1} = 0; % do not orthogonalise them  
           
           [RPE,psi] = get_latents(subj_original_idx, run);

           multi.pmod(1).name{1} = 'psi';
           multi.pmod(1).param{1} = psi';
           multi.pmod(1).poly{1} = 1; 

            % nuisance @ reaction onset
            %
            multi.names{2} = 'reaction_onset';
            multi.onsets{2} = reaction_onsets;
            multi.durations{2} = zeros(size(multi.onsets{2}));

            % nuisance @ trial onset
            %
            multi.names{3} = 'trial_onset';
            multi.onsets{3} = trial_onsets;
            multi.durations{3} = zeros(size(multi.onsets{3}));
            
            % nuisance @ feedback onset
            %
            multi.names{4} = 'feedback_onset';
            multi.onsets{4} = feedback_onsets;
            multi.durations{4} = zeros(size(multi.onsets{4}));
                        
               
            if sum(timeouts_latent_guess) > 0
               multi.names{5} = 'feedback_onset_timeouts';
               multi.onsets{5} = feedback_onsets(timeouts_latent_guess); % timeouts only
               multi.durations{5} = zeros(size(multi.onsets{5}));
            end
            
                   
      % RPE low and high
        % @ feedback
        % nuisance @ 
    
       
        case 19 
            
            event_index = 1;
            
            [RPE,psi] = get_latents(subj_original_idx, run);
            
            psi_low = psi < 0.52;
            psi_high = ~psi_low;
            
            psi_low_1 = logical(zeros(size(feedback_onsets)));
            psi_low_1(~timeouts_latent_guess) = psi_low;
            
            psi_high_1 = logical(zeros(size(feedback_onsets)));
            psi_high_1(~timeouts_latent_guess) = psi_high;
            
                    
            save debug.mat
            
            assert(sum(psi_low_1) + sum(psi_high_1) <= length(feedback_onsets))
            assert(length(psi_low_1) == length(feedback_onsets))
            assert(length(psi_high_1) == length(feedback_onsets))
           
             if any(psi_low_1)         
                % psi low
                multi.names{event_index} = ['psi_low'];
                multi.onsets{event_index} = feedback_onsets(psi_low_1);
                multi.durations{event_index} = zeros(size(multi.onsets{event_index})); %impulse regressor (as opposed to boxcar regressor)
                
                multi.pmod(event_index).name{1} = 'RPE';
                multi.pmod(event_index).param{1} = RPE(psi_low)';
                multi.pmod(event_index).poly{1} = 1; 
                
                event_index = event_index +1;
             end   
            
             if any(psi_high_1)         
                % psi high
                multi.names{event_index} = ['psi_high'];
                multi.onsets{event_index} = feedback_onsets(psi_high_1);
                multi.durations{event_index} = zeros(size(multi.onsets{event_index})); %impulse regressor (as opposed to boxcar regressor)
                
                multi.pmod(event_index).name{1} = 'RPE';
                multi.pmod(event_index).param{1} = RPE(psi_high)';
                multi.pmod(event_index).poly{1} = 1; 
                
                event_index = event_index +1;
             end
             
                 
            % nuisance @ reaction onset
            %
            multi.names{event_index} = 'reaction_onset';
            multi.onsets{event_index} = reaction_onsets;
            multi.durations{event_index} = zeros(size(multi.onsets{event_index}));
            event_index = event_index +1;

            % nuisance @ trial onset
            %
            multi.names{event_index} = 'trial_onset';
            multi.onsets{event_index} = trial_onsets;
            multi.durations{event_index} = zeros(size(multi.onsets{event_index}));
            event_index = event_index +1;

           % nuisance @ latent onset
            %
            multi.names{event_index} = 'latent_onset';
            multi.onsets{event_index} = latent_onsets;
            multi.durations{event_index} = zeros(size(multi.onsets{event_index}));
            
            % nuisance @ latent offset
            %
            multi.names{event_index} = 'latent_offset';
            multi.onsets{event_index} = latent_offsets;
            multi.durations{event_index} = zeros(size(multi.onsets{event_index}));
         
               
            if sum(timeouts) > 0
               multi.names{event_index} = 'feedback_onset_timeouts';
               multi.onsets{event_index} = feedback_onsets(timeouts); % timeouts only
               multi.durations{event_index} = zeros(size(multi.onsets{event_index}));
            end
                  
            
         % RPE * 1-PSI (multiplicative) @ feedback onset
        % nuisance @ trial and choice onset 
        % 
        case 20 
            % 
             
    
            multi.names{1} = condition;
            multi.onsets{1} = feedback_onsets(~timeouts_latent_guess);
            multi.durations{1} = zeros(size(multi.onsets{1})); %impulse regressor (as opposed to boxcar regressor)

           multi.orth{1} = 0; % do not orthogonalise them  
           
           [RPE,psi] = get_latents(subj_original_idx, run);
           
           multi.pmod(1).name{1} = 'RPE_1_minus_psi';
           multi.pmod(1).param{1} = RPE' .* (1-psi)';
           multi.pmod(1).poly{1} = 1;   
           

            % nuisance @ reaction onset
            %
            multi.names{2} = 'reaction_onset';
            multi.onsets{2} = reaction_onsets;
            multi.durations{2} = zeros(size(multi.onsets{2}));

            % nuisance @ trial onset
            %
            multi.names{3} = 'trial_onset';
            multi.onsets{3} = trial_onsets;
            multi.durations{3} = zeros(size(multi.onsets{3}));
            
            % nuisance @ latent onset
            %
            multi.names{4} = 'latent_onset';
            multi.onsets{4} = latent_onsets;
            multi.durations{4} = zeros(size(multi.onsets{4}));
            
            % nuisance @ latent offset
            %
            multi.names{5} = 'latent_offset';
            multi.onsets{5} = latent_offsets;
            multi.durations{5} = zeros(size(multi.onsets{5}));
               
            if sum(timeouts_latent_guess) > 0
               multi.names{6} = 'feedback_onset_timeouts';
               multi.onsets{6} = feedback_onsets(timeouts_latent_guess); % timeouts only
               multi.durations{6} = zeros(size(multi.onsets{6}));
            end     
            
            
            % RPE * (1-PSI) (pavlovian) and RPE*psi (instrumental) @ feedback onset
        % nuisance @ trial and choice onset 
        % 
        case 21 
            % 
            %
            
    
            multi.names{1} = condition;
            multi.onsets{1} = feedback_onsets(~timeouts_latent_guess);
            multi.durations{1} = zeros(size(multi.onsets{1})); %impulse regressor (as opposed to boxcar regressor)

           multi.orth{1} = 0; % do not orthogonalise them  
           
           [RPE,psi] = get_latents(subj_original_idx, run);
           
           multi.pmod(1).name{1} = 'RPEpsi';
           multi.pmod(1).param{1} = RPE' .* psi';
           multi.pmod(1).poly{1} = 1; 
           
           multi.pmod(1).name{2} = 'RPE_1_minus_psi';
           multi.pmod(1).param{2} = RPE' .* (1-psi)';
           multi.pmod(1).poly{2} = 1;   
           

            % nuisance @ reaction onset
            %
            multi.names{2} = 'reaction_onset';
            multi.onsets{2} = reaction_onsets;
            multi.durations{2} = zeros(size(multi.onsets{2}));

            % nuisance @ feedback onset
            %
            multi.names{3} = 'trial_onset';
            multi.onsets{3} = trial_onsets;
            multi.durations{3} = zeros(size(multi.onsets{3}));
            
            % nuisance @ latent onset
            %
            multi.names{4} = 'latent_onset';
            multi.onsets{4} = latent_onsets;
            multi.durations{4} = zeros(size(multi.onsets{4}));
            
            % nuisance @ latent offset
            %
            multi.names{5} = 'latent_offset';
            multi.onsets{5} = latent_offsets;
            multi.durations{5} = zeros(size(multi.onsets{5}));
            
               
            if sum(timeouts_latent_guess) > 0
               multi.names{6} = 'feedback_onset_timeouts';
               multi.onsets{6} = feedback_onsets(timeouts_latent_guess); % timeouts only
               multi.durations{6} = zeros(size(multi.onsets{6}));
            end     
                    
               
            
 
         % Psi only @ latent onset
      
        % 
        case 22
            % 
            %
            
    
            multi.names{1} = condition;
            multi.onsets{1} = latent_onsets(~timeouts_latent_guess);
            multi.durations{1} = zeros(size(multi.onsets{1})); %impulse regressor (as opposed to boxcar regressor)

           multi.orth{1} = 0; % do not orthogonalise them  
           
           [RPE,psi] = get_latents(subj_original_idx, run);

           multi.pmod(1).name{1} = 'psi';
           multi.pmod(1).param{1} = psi';
           multi.pmod(1).poly{1} = 1; 

            % nuisance @ reaction onset
            %
            multi.names{2} = 'reaction_onset';
            multi.onsets{2} = reaction_onsets;
            multi.durations{2} = zeros(size(multi.onsets{2}));

            % nuisance @ feedback onset
            %
            multi.names{3} = 'trial_onset';
            multi.onsets{3} = trial_onsets;
            multi.durations{3} = zeros(size(multi.onsets{3}));
            
            % nuisance @ feedback onset
            %
            multi.names{4} = 'feedback_onset';
            multi.onsets{4} = feedback_onsets;
            multi.durations{4} = zeros(size(multi.onsets{4}));
            
            % nuisance @ latent offset
            %
            multi.names{5} = 'latent_offset';
            multi.onsets{5} = latent_offsets;
            multi.durations{5} = zeros(size(multi.onsets{5}));
            
               
            if sum(timeouts_latent_guess) > 0
               multi.names{6} = 'feedback_onset_timeouts';
               multi.onsets{6} = feedback_onsets(timeouts_latent_guess); % timeouts only
               multi.durations{6} = zeros(size(multi.onsets{6}));
            end
                       
            
% yes > no at subjective latent guess
           case 23 
            
            event_index = 1;
            
            save debug.mat; 
            
            assert(sum(yes) + sum(no) <= length(latent_onsets))
            assert(length(yes) == length(latent_onsets))
            assert(length(no) == length(latent_onsets))
            
             if any(yes)         
                % latent guess = yes
                multi.names{event_index} = ['yes'];
                multi.onsets{event_index} = latent_onsets(yes);
                multi.durations{event_index} = zeros(size(multi.onsets{event_index})); %impulse regressor (as opposed to boxcar regressor)
                event_index = event_index +1;
             end
             
             if any(no)         
                % latent guess = no
                multi.names{event_index} = ['no'];
                multi.onsets{event_index} = latent_onsets(no);
                multi.durations{event_index} = zeros(size(multi.onsets{event_index})); %impulse regressor (as opposed to boxcar regressor)
                event_index = event_index +1;
             end
             
             % nuisance @ reaction onset
            %
            multi.names{event_index} = 'reaction_onset';
            multi.onsets{event_index} = reaction_onsets;
            multi.durations{event_index} = zeros(size(multi.onsets{event_index}));
            event_index = event_index +1;

            % nuisance @ trial onset
            %
            multi.names{event_index} = 'trial_onset';
            multi.onsets{event_index} = trial_onsets;
            multi.durations{event_index} = zeros(size(multi.onsets{event_index}));
            event_index = event_index +1;
            
            % nuisance @ fb onset
            %
            multi.names{event_index} = 'feedback_onset';
            multi.onsets{event_index} = feedback_onsets;
            multi.durations{event_index} = zeros(size(multi.onsets{event_index}));
            
            % nuisance @ latent offset
            %
            multi.names{event_index} = 'latent_offset';
            multi.onsets{event_index} = latent_offsets;
            multi.durations{event_index} = zeros(size(multi.onsets{event_index}));
         
               
            if sum(timeouts) > 0
               multi.names{event_index} = 'feedback_onset_timeouts';
               multi.onsets{event_index} = feedback_onsets(timeouts); % timeouts only
               multi.durations{event_index} = zeros(size(multi.onsets{event_index}));
            end         
               
            
        % Functional connectivity - from psi GLM seed region (GLM 11)
        case 24    
            
            multi.names{1} = condition;
            multi.onsets{1} = feedback_onsets(~timeouts_latent_guess);
            multi.durations{1} = zeros(size(multi.onsets{1})); %impulse regressor (as opposed to boxcar regressor)

           multi.orth{1} = 0; % do not orthogonalise them  
           
            %pg 133 from Poldrack book for beta series correlation 
           
           [RPE,psi] = get_latents(subj_original_idx, run);
           %disp(psi)
           
           MidTemp = getBetaSeries(optCon_expt(),7,subj, run, 'feedback_onset', 'masks/glm11_psi_cluster_t=6.282_extent=220_roi=Temporal_Mid_R_peak=[54_-26_-12].nii');
           MidTemp = MidTemp(~timeouts_latent_guess);
%         
%            multi.pmod(1).name{1} = 'psi'; % this is the X
%            multi.pmod(1).param{1} = psi';
%            multi.pmod(1).poly{1} = 1;
           
           multi.pmod(1).name{1} = 'MidTemp'; %this is the R
           multi.pmod(1).param{1} = MidTemp';
           multi.pmod(1).poly{1} = 1;
           

            % nuisance @ reaction onset
            %
            multi.names{2} = 'reaction_onset';
            multi.onsets{2} = reaction_onsets;
            multi.durations{2} = zeros(size(multi.onsets{2}));

            % nuisance @ trial onset
            %
            multi.names{3} = 'trial_onset';
            multi.onsets{3} = trial_onsets;
            multi.durations{3} = zeros(size(multi.onsets{3}));
            
            % nuisance @ latent onset
            %
            multi.names{4} = 'latent_onset';
            multi.onsets{4} = latent_onsets;
            multi.durations{4} = zeros(size(multi.onsets{4}));
            
            % nuisance @ latent offset
            %
            multi.names{5} = 'latent_offset';
            multi.onsets{5} = latent_offsets;
            multi.durations{5} = zeros(size(multi.onsets{5}));
            
               
            if sum(timeouts_latent_guess) > 0
               multi.names{6} = 'feedback_onset_timeouts';
               multi.onsets{6} = feedback_onsets(timeouts_latent_guess); % timeouts only
               multi.durations{6} = zeros(size(multi.onsets{6}));
            end     
            
            
            
            
             
        % RPE * PSI (multiplicative) and RPE and PSI @ feedback onset
        % nuisance @ trial and choice onset 
        % 
        case 25 
            % 
            %
            
    
            multi.names{1} = condition;
            multi.onsets{1} = feedback_onsets(~timeouts_latent_guess);
            multi.durations{1} = zeros(size(multi.onsets{1})); %impulse regressor (as opposed to boxcar regressor)

           multi.orth{1} = 0; % do not orthogonalise them  
           
           [RPE,psi] = get_latents(subj_original_idx, run);
           
           multi.pmod(1).name{1} = 'RPEpsi';
           multi.pmod(1).param{1} = RPE' .* psi';
           multi.pmod(1).poly{1} = 1; 
           
           multi.pmod(1).name{2} = 'RPE';
           multi.pmod(1).param{2} = RPE';
           multi.pmod(1).poly{2} = 1;   
           
           multi.pmod(1).name{3} = 'psi';
           multi.pmod(1).param{3} = psi';
           multi.pmod(1).poly{3} = 1;  
           

            % nuisance @ reaction onset
            %
            multi.names{2} = 'reaction_onset';
            multi.onsets{2} = reaction_onsets;
            multi.durations{2} = zeros(size(multi.onsets{2}));

            % nuisance @ feedback onset
            %
            multi.names{3} = 'trial_onset';
            multi.onsets{3} = trial_onsets;
            multi.durations{3} = zeros(size(multi.onsets{3}));
            
            % nuisance @ latent onset
            %
            multi.names{4} = 'latent_onset';
            multi.onsets{4} = latent_onsets;
            multi.durations{4} = zeros(size(multi.onsets{4}));
            
            % nuisance @ latent offset
            %
            multi.names{5} = 'latent_offset';
            multi.onsets{5} = latent_offsets;
            multi.durations{5} = zeros(size(multi.onsets{5}));
            
               
            if sum(timeouts_latent_guess) > 0
               multi.names{6} = 'feedback_onset_timeouts';
               multi.onsets{6} = feedback_onsets(timeouts_latent_guess); % timeouts only
               multi.durations{6} = zeros(size(multi.onsets{6}));
            end     
            
       % yes > no at subjective latent guess offset (button press)
           case 26 
            
            event_index = 1;
            
            assert(sum(yes) + sum(no) <= length(latent_offsets))
            assert(length(yes) == length(latent_offsets))
            assert(length(no) == length(latent_offsets))
            
             if any(yes)         
                % latent guess = yes
                multi.names{event_index} = ['yes'];
                multi.onsets{event_index} = latent_offsets(yes);
                multi.durations{event_index} = zeros(size(multi.onsets{event_index})); %impulse regressor (as opposed to boxcar regressor)
                event_index = event_index +1;
             end
             
             if any(no)         
                % latent guess = no
                multi.names{event_index} = ['no'];
                multi.onsets{event_index} = latent_offsets(no);
                multi.durations{event_index} = zeros(size(multi.onsets{event_index})); %impulse regressor (as opposed to boxcar regressor)
                event_index = event_index +1;
             end
             
             % nuisance @ reaction onset
            %
            multi.names{event_index} = 'reaction_onset';
            multi.onsets{event_index} = reaction_onsets;
            multi.durations{event_index} = zeros(size(multi.onsets{event_index}));
            event_index = event_index +1;

            % nuisance @ trial onset
            %
            multi.names{event_index} = 'trial_onset';
            multi.onsets{event_index} = trial_onsets;
            multi.durations{event_index} = zeros(size(multi.onsets{event_index}));
            event_index = event_index +1;
            
            % nuisance @ fb onset
            %
            multi.names{event_index} = 'feedback_onset';
            multi.onsets{event_index} = feedback_onsets;
            multi.durations{event_index} = zeros(size(multi.onsets{event_index}));
            
            % nuisance @ latent onset
            %
            multi.names{event_index} = 'latent_onset';
            multi.onsets{event_index} = latent_onsets;
            multi.durations{event_index} = zeros(size(multi.onsets{event_index}));
         
               
            if sum(timeouts) > 0
               multi.names{event_index} = 'feedback_onset_timeouts';
               multi.onsets{event_index} = feedback_onsets(timeouts); % timeouts only
               multi.durations{event_index} = zeros(size(multi.onsets{event_index}));
            end    
            
            
      % Yes vs. No in Benev vs. Adv
        % Latent guess @ latent onset 
        % nuisance @ reaction and trial onset and latent
        % offset(button press)

        
        case 27 
            
            event_index = 1;
            
  
            if any(yes)         
                % yeses for condition
                multi.names{event_index} = ['yes_',condition];
                multi.onsets{event_index} = latent_onsets(yes);
                multi.durations{event_index} = zeros(size(multi.onsets{event_index})); %impulse regressor (as opposed to boxcar regressor)
                event_index = event_index +1;
             end
           
            if any(no)
                % nos for condition
                multi.names{event_index} = ['no_',condition];
                multi.onsets{event_index} = latent_onsets(no);
                multi.durations{event_index} = zeros(size(multi.onsets{event_index})); %impulse regressor (as opposed to boxcar regressor)
                event_index = event_index +1;
            end
             

            % nuisance @ reaction onset
            %
            multi.names{event_index} = 'reaction_onset';
            multi.onsets{event_index} = reaction_onsets;
            multi.durations{event_index} = zeros(size(multi.onsets{event_index}));
            event_index = event_index +1;

            % nuisance @ trial onset
            %
            multi.names{event_index} = 'trial_onset';
            multi.onsets{event_index} = trial_onsets;
            multi.durations{event_index} = zeros(size(multi.onsets{event_index}));
            event_index = event_index +1;
            
            % nuisance @ latent onset
            %
            multi.names{event_index} = 'feedback_onset';
            multi.onsets{event_index} = feedback_onsets;
            multi.durations{event_index} = zeros(size(multi.onsets{event_index}));
            
            % nuisance @ latent offset
            %
            multi.names{event_index} = 'latent_offset';
            multi.onsets{event_index} = latent_offsets;
            multi.durations{event_index} = zeros(size(multi.onsets{event_index}));
         
               
            if sum(timeouts) > 0
               multi.names{event_index} = 'latent_onset_timeouts';
               multi.onsets{event_index} = latent_onsets(timeouts); % timeouts only
               multi.durations{event_index} = zeros(size(multi.onsets{event_index}));
            end
            
            
       % delta and beta @ feedback onset
        % nuisance @ trial and choice onset 
        % 
        case 28
            % 
            %
            event_index = 1;
                        
    
            multi.names{1} = condition;
            multi.onsets{1} = feedback_onsets(~timeouts_latent_guess);
            multi.durations{1} = zeros(size(multi.onsets{1})); %impulse regressor (as opposed to boxcar regressor)

           multi.orth{1} = 0; % do not orthogonalise them  
           
           delta = wins - losses;
           delta = delta(~timeouts_latent_guess);
           
           b = double(no(~timeouts_latent_guess));
           
           
           multi.pmod(1).name{1} = 'delta';
           multi.pmod(1).param{1} = delta;
           multi.pmod(1).poly{1} = 1;  
           
           multi.pmod(1).name{2} = 'b';
           multi.pmod(1).param{2} = b;
           multi.pmod(1).poly{2} = 1; 
           
           multi.pmod(1).name{3} = 'delta_b';
           multi.pmod(1).param{3} = delta .* b;
           multi.pmod(1).poly{3} = 1; 
           

            % nuisance @ reaction onset
            %
            multi.names{2} = 'reaction_onset';
            multi.onsets{2} = reaction_onsets;
            multi.durations{2} = zeros(size(multi.onsets{2}));

            % nuisance @ trial onset
            %
            multi.names{3} = 'trial_onset';
            multi.onsets{3} = trial_onsets;
            multi.durations{3} = zeros(size(multi.onsets{3}));
            
            % nuisance @ latent onset
            %
            multi.names{4} = 'latent_onset';
            multi.onsets{4} = latent_onsets;
            multi.durations{4} = zeros(size(multi.onsets{4}));
            
            % nuisance @ latent offset
            %
            multi.names{5} = 'latent_offset';
            multi.onsets{5} = latent_offsets;
            multi.durations{5} = zeros(size(multi.onsets{5}));

            
               
            if sum(timeouts_latent_guess) > 0
               multi.names{6} = 'feedback_onset_timeouts';
               multi.onsets{6} = feedback_onsets(timeouts_latent_guess); % timeouts only
               multi.durations{6} = zeros(size(multi.onsets{6}));
            end     
            
            
           
            
             % @feedback
        case 29    
                       
            event_index = 1;
                      
            if any(~timeouts_latent_guess)
                %all trials
                multi.names{event_index} = 'all_trials';
                multi.onsets{event_index} = feedback_onsets(~timeouts_latent_guess);
                multi.durations{event_index} = zeros(size(multi.onsets{event_index})); %impulse regressor (as opposed to boxcar regressor)
                event_index = event_index +1;
            end
%             
            
            if any(wins & ~timeouts_latent_guess)
            %win trials
                multi.names{event_index} = 'wins';
                multi.onsets{event_index} = feedback_onsets(wins & ~timeouts_latent_guess);
                multi.durations{event_index} = zeros(size(multi.onsets{event_index})); %impulse regressor (as opposed to boxcar regressor)
                event_index = event_index +1;
            end
%             
            if any(adversarial & ~timeouts_latent_guess)
            %adversarial trials
                multi.names{event_index} = 'adversarial';
                multi.onsets{event_index} = feedback_onsets(adversarial & ~timeouts_latent_guess);
                multi.durations{event_index} = zeros(size(multi.onsets{event_index})); %impulse regressor (as opposed to boxcar regressor)
                event_index = event_index +1;
            end
            
            if any(wins & adversarial & ~timeouts_latent_guess)
                %win adversarial trials
                multi.names{event_index} = 'wins_adversarial';
                multi.onsets{event_index} = feedback_onsets(wins & adversarial & ~timeouts_latent_guess);
                multi.durations{event_index} = zeros(size(multi.onsets{event_index})); %impulse regressor (as opposed to boxcar regressor)
                event_index = event_index +1;
            end
           

           multi.orth{1} = 0; % do not orthogonalise them            
                
            
            % nuisance @ reaction onset
            %
            multi.names{event_index} = 'reaction_onset';
            multi.onsets{event_index} = reaction_onsets;
            multi.durations{event_index} = zeros(size(multi.onsets{event_index}));
            event_index = event_index +1;
         

            % nuisance @ trial onset
            %
            multi.names{event_index} = 'trial_onset';
            multi.onsets{event_index} = trial_onsets;
            multi.durations{event_index} = zeros(size(multi.onsets{event_index}));
            event_index = event_index +1;
       

           % nuisance @ latent onset
            %
            multi.names{event_index} = 'latent_onset';
            multi.onsets{event_index} = latent_onsets;
            multi.durations{event_index} = zeros(size(multi.onsets{event_index}));
            event_index = event_index +1;
            
            % nuisance @ latent offset
            %
            multi.names{event_index} = 'latent_offset';
            multi.onsets{event_index} = latent_offsets;
            multi.durations{event_index} = zeros(size(multi.onsets{event_index}));
            event_index = event_index +1;
        
            
            if sum(timeouts_latent_guess) > 0
               multi.names{event_index} = 'feedback_onset_timeouts';
               multi.onsets{event_index} = feedback_onsets(timeouts_latent_guess); % timeouts only
               multi.durations{event_index} = zeros(size(multi.onsets{event_index}));
               event_index = event_index +1;
            end   
            
               
            
            case 30    
                       
            event_index = 1;
            
        
            
%             if any(feedback_onsets)
%                 %all trials
%                 multi.names{event_index} = 'all_trials';
%                 multi.onsets{event_index} = feedback_onsets(~timeouts_latent_guess);
%                 multi.durations{event_index} = zeros(size(multi.onsets{event_index})); %impulse regressor (as opposed to boxcar regressor)
%                 event_index = event_index +1;
%             end
%             
            
            if any(wins)
            %win trials
                multi.names{event_index} = 'wins';
                multi.onsets{event_index} = feedback_onsets(wins);
                multi.durations{event_index} = zeros(size(multi.onsets{event_index})); %impulse regressor (as opposed to boxcar regressor)
                event_index = event_index +1;
            end
%             
            if any(adversarial)
            %adversarial trials
                multi.names{event_index} = 'adversarial';
                multi.onsets{event_index} = feedback_onsets(adversarial);
                multi.durations{event_index} = zeros(size(multi.onsets{event_index})); %impulse regressor (as opposed to boxcar regressor)
                event_index = event_index +1;
            end
            
            if any(wins_adversarial)
                %win adversarial trials
                multi.names{event_index} = 'wins_adversarial';
                multi.onsets{event_index} = feedback_onsets(wins_adversarial);
                multi.durations{event_index} = zeros(size(multi.onsets{event_index})); %impulse regressor (as opposed to boxcar regressor)
                event_index = event_index +1;
            end
           

           multi.orth{1} = 0; % do not orthogonalise them            
                
            
            % nuisance @ reaction onset
            %
            multi.names{event_index} = 'reaction_onset';
            multi.onsets{event_index} = reaction_onsets;
            multi.durations{event_index} = zeros(size(multi.onsets{event_index}));
            event_index = event_index +1;
         

            % nuisance @ trial onset
            %
            multi.names{event_index} = 'trial_onset';
            multi.onsets{event_index} = trial_onsets;
            multi.durations{event_index} = zeros(size(multi.onsets{event_index}));
            event_index = event_index +1;
       

           % nuisance @ latent onset
            %
            multi.names{event_index} = 'latent_onset';
            multi.onsets{event_index} = latent_onsets;
            multi.durations{event_index} = zeros(size(multi.onsets{event_index}));
            event_index = event_index +1;
            
            % nuisance @ latent offset
            %
            multi.names{event_index} = 'latent_offset';
            multi.onsets{event_index} = latent_offsets;
            multi.durations{event_index} = zeros(size(multi.onsets{event_index}));
            event_index = event_index +1;
        
            
            if sum(timeouts_latent_guess) > 0
               multi.names{event_index} = 'feedback_onset_timeouts';
               multi.onsets{event_index} = feedback_onsets(timeouts_latent_guess); % timeouts only
               multi.durations{event_index} = zeros(size(multi.onsets{event_index}));
               event_index = event_index +1;
            end   
            
               
            
            case 31    
                       
            event_index = 1;
            
        
            
%             if any(feedback_onsets)
%                 %all trials
%                 multi.names{event_index} = 'all_trials';
%                 multi.onsets{event_index} = feedback_onsets(~timeouts_latent_guess);
%                 multi.durations{event_index} = zeros(size(multi.onsets{event_index})); %impulse regressor (as opposed to boxcar regressor)
%                 event_index = event_index +1;
%             end
%             
            
%             if any(wins)
%             %win trials
%                 multi.names{event_index} = 'wins';
%                 multi.onsets{event_index} = feedback_onsets(wins);
%                 multi.durations{event_index} = zeros(size(multi.onsets{event_index})); %impulse regressor (as opposed to boxcar regressor)
%                 event_index = event_index +1;
%             end
%             
%             if any(adversarial)
%             %adversarial trials
%                 multi.names{event_index} = 'adversarial';
%                 multi.onsets{event_index} = feedback_onsets(adversarial);
%                 multi.durations{event_index} = zeros(size(multi.onsets{event_index})); %impulse regressor (as opposed to boxcar regressor)
%                 event_index = event_index +1;
%             end
            
            if any(wins_adversarial)
                %win adversarial trials
                multi.names{event_index} = 'wins_adversarial';
                multi.onsets{event_index} = feedback_onsets(wins_adversarial);
                multi.durations{event_index} = zeros(size(multi.onsets{event_index})); %impulse regressor (as opposed to boxcar regressor)
                event_index = event_index +1;
            end
           

           multi.orth{1} = 0; % do not orthogonalise them            
                
            
            % nuisance @ reaction onset
            %
            multi.names{event_index} = 'reaction_onset';
            multi.onsets{event_index} = reaction_onsets;
            multi.durations{event_index} = zeros(size(multi.onsets{event_index}));
            event_index = event_index +1;
         

            % nuisance @ trial onset
            %
            multi.names{event_index} = 'trial_onset';
            multi.onsets{event_index} = trial_onsets;
            multi.durations{event_index} = zeros(size(multi.onsets{event_index}));
            event_index = event_index +1;
       

           % nuisance @ latent onset
            %
            multi.names{event_index} = 'latent_onset';
            multi.onsets{event_index} = latent_onsets;
            multi.durations{event_index} = zeros(size(multi.onsets{event_index}));
            event_index = event_index +1;
            
            % nuisance @ latent offset
            %
            multi.names{event_index} = 'latent_offset';
            multi.onsets{event_index} = latent_offsets;
            multi.durations{event_index} = zeros(size(multi.onsets{event_index}));
            event_index = event_index +1;
        
            
            if sum(timeouts_latent_guess) > 0
               multi.names{event_index} = 'feedback_onset_timeouts';
               multi.onsets{event_index} = feedback_onsets(timeouts_latent_guess); % timeouts only
               multi.durations{event_index} = zeros(size(multi.onsets{event_index}));
               event_index = event_index +1;
            end
            
            
            
            
            
        case 32
                       
            event_index = 1;
            
        
            
%             if any(feedback_onsets)
%                 %all trials
%                 multi.names{event_index} = 'all_trials';
%                 multi.onsets{event_index} = feedback_onsets(~timeouts_latent_guess);
%                 multi.durations{event_index} = zeros(size(multi.onsets{event_index})); %impulse regressor (as opposed to boxcar regressor)
%                 event_index = event_index +1;
%             end
%             
            
            if any(wins)
            %win trials
                multi.names{event_index} = 'wins';
                multi.onsets{event_index} = feedback_onsets(wins);
                multi.durations{event_index} = zeros(size(multi.onsets{event_index})); %impulse regressor (as opposed to boxcar regressor)
                event_index = event_index +1;
            end
%             
            if any(adversarial)
            %adversarial trials
                multi.names{event_index} = 'adversarial';
                multi.onsets{event_index} = feedback_onsets(adversarial);
                multi.durations{event_index} = zeros(size(multi.onsets{event_index})); %impulse regressor (as opposed to boxcar regressor)
                event_index = event_index +1;
            end
            
            if any(wins_adversarial)
                %win adversarial trials
                multi.names{event_index} = 'wins_adversarial';
                multi.onsets{event_index} = feedback_onsets(wins_adversarial);
                multi.durations{event_index} = zeros(size(multi.onsets{event_index})); %impulse regressor (as opposed to boxcar regressor)
                event_index = event_index +1;
            end
           

           multi.orth{1} = 0; % do not orthogonalise them            
                
            
            % nuisance @ reaction onset
            %
            multi.names{event_index} = 'reaction_onset';
            multi.onsets{event_index} = reaction_onsets;
            multi.durations{event_index} = zeros(size(multi.onsets{event_index}));
            event_index = event_index +1;
         

            % nuisance @ trial onset
            %
            multi.names{event_index} = 'trial_onset';
            multi.onsets{event_index} = trial_onsets;
            multi.durations{event_index} = zeros(size(multi.onsets{event_index}));
            event_index = event_index +1;
       

           % nuisance @ latent onset
            %
            multi.names{event_index} = 'latent_onset';
            multi.onsets{event_index} = latent_onsets;
            multi.durations{event_index} = zeros(size(multi.onsets{event_index}));
            event_index = event_index +1;
            
            % nuisance @ latent offset
            %
            multi.names{event_index} = 'latent_offset';
            multi.onsets{event_index} = latent_offsets;
            multi.durations{event_index} = zeros(size(multi.onsets{event_index}));
            event_index = event_index +1;
        
            
            if sum(timeouts_latent_guess) > 0
               multi.names{event_index} = 'feedback_onset_timeouts';
               multi.onsets{event_index} = feedback_onsets(timeouts_latent_guess); % timeouts only
               multi.durations{event_index} = zeros(size(multi.onsets{event_index}));
               event_index = event_index +1;
            end   
                      
              % @feedback
        case 33    
                       
            event_index = 1;
            
        
%             
            if any(feedback_onsets)
                %all trials
                multi.names{event_index} = 'all_trials';
                multi.onsets{event_index} = feedback_onsets(~timeouts_latent_guess);
                multi.durations{event_index} = zeros(size(multi.onsets{event_index})); %impulse regressor (as opposed to boxcar regressor)
                event_index = event_index +1;
            end
%             
            
            if any(wins)
            %win trials
                multi.names{event_index} = 'wins';
                multi.onsets{event_index} = feedback_onsets(wins);
                multi.durations{event_index} = zeros(size(multi.onsets{event_index})); %impulse regressor (as opposed to boxcar regressor)
                event_index = event_index +1;
            end
%             
            if any(adversarial)
            %adversarial trials
                multi.names{event_index} = 'adversarial';
                multi.onsets{event_index} = feedback_onsets(adversarial);
                multi.durations{event_index} = zeros(size(multi.onsets{event_index})); %impulse regressor (as opposed to boxcar regressor)
                event_index = event_index +1;
            end
            
%             if any(wins_adversarial)
%                 %win adversarial trials
%                 multi.names{event_index} = 'wins_adversarial';
%                 multi.onsets{event_index} = feedback_onsets(wins_adversarial);
%                 multi.durations{event_index} = zeros(size(multi.onsets{event_index})); %impulse regressor (as opposed to boxcar regressor)
%                 event_index = event_index +1;
%             end
           

           multi.orth{1} = 0; % do not orthogonalise them            
                
            
            % nuisance @ reaction onset
            %
            multi.names{event_index} = 'reaction_onset';
            multi.onsets{event_index} = reaction_onsets;
            multi.durations{event_index} = zeros(size(multi.onsets{event_index}));
            event_index = event_index +1;
         

            % nuisance @ trial onset
            %
            multi.names{event_index} = 'trial_onset';
            multi.onsets{event_index} = trial_onsets;
            multi.durations{event_index} = zeros(size(multi.onsets{event_index}));
            event_index = event_index +1;
       

           % nuisance @ latent onset
            %
            multi.names{event_index} = 'latent_onset';
            multi.onsets{event_index} = latent_onsets;
            multi.durations{event_index} = zeros(size(multi.onsets{event_index}));
            event_index = event_index +1;
            
            % nuisance @ latent offset
            %
            multi.names{event_index} = 'latent_offset';
            multi.onsets{event_index} = latent_offsets;
            multi.durations{event_index} = zeros(size(multi.onsets{event_index}));
            event_index = event_index +1;
        
            
            if sum(timeouts_latent_guess) > 0
               multi.names{event_index} = 'feedback_onset_timeouts';
               multi.onsets{event_index} = feedback_onsets(timeouts_latent_guess); % timeouts only
               multi.durations{event_index} = zeros(size(multi.onsets{event_index}));
               event_index = event_index +1;
            end   
            
            
            
   % @feedback --- original way
        case 34    
                       
            event_index = 1;
            
        
            
            if any(feedback_onsets & ~timeouts_latent_guess)
                %all trials
                multi.names{event_index} = 'all_trials';
                multi.onsets{event_index} = feedback_onsets(~timeouts_latent_guess);
                multi.durations{event_index} = zeros(size(multi.onsets{event_index})); %impulse regressor (as opposed to boxcar regressor)
                event_index = event_index +1;
            end
            
            
            if any(wins & ~timeouts_latent_guess)
            %win trials
                multi.names{event_index} = 'wins';
                multi.onsets{event_index} = feedback_onsets(wins & ~timeouts_latent_guess);
                multi.durations{event_index} = zeros(size(multi.onsets{event_index})); %impulse regressor (as opposed to boxcar regressor)
                event_index = event_index +1;
            end
            
            if any(adversarial & ~timeouts_latent_guess)
            %adversarial trials
                multi.names{event_index} = 'adversarial';
                multi.onsets{event_index} = feedback_onsets(adversarial & ~timeouts_latent_guess);
                multi.durations{event_index} = zeros(size(multi.onsets{event_index})); %impulse regressor (as opposed to boxcar regressor)
                event_index = event_index +1;
            end
            
            if any(wins & adversarial & ~timeouts_latent_guess)
                %win adversarial trials
                multi.names{event_index} = 'wins_adversarial';
                multi.onsets{event_index} = feedback_onsets(wins & adversarial & ~timeouts_latent_guess);
                multi.durations{event_index} = zeros(size(multi.onsets{event_index})); %impulse regressor (as opposed to boxcar regressor)
                event_index = event_index +1;
            end
           

           multi.orth{1} = 0; % do not orthogonalise them            
                
            
            %nuisance @ reaction onset
            
            multi.names{event_index} = 'reaction_onset';
            multi.onsets{event_index} = reaction_onsets;
            multi.durations{event_index} = zeros(size(multi.onsets{event_index}));
            event_index = event_index +1;
         

            %nuisance @ trial onset
            
            multi.names{event_index} = 'trial_onset';
            multi.onsets{event_index} = trial_onsets;
            multi.durations{event_index} = zeros(size(multi.onsets{event_index}));
            event_index = event_index +1;
       

           %nuisance @ latent onset
            
            multi.names{event_index} = 'latent_onset';
            multi.onsets{event_index} = latent_onsets;
            multi.durations{event_index} = zeros(size(multi.onsets{event_index}));
            event_index = event_index +1;
            
           % nuisance @ latent offset
            
            multi.names{event_index} = 'latent_offset';
            multi.onsets{event_index} = latent_offsets;
            multi.durations{event_index} = zeros(size(multi.onsets{event_index}));
            event_index = event_index +1;
        
            
            if sum(timeouts_latent_guess) > 0
               multi.names{event_index} = 'feedback_onset_timeouts';
               multi.onsets{event_index} = feedback_onsets(timeouts_latent_guess); % timeouts only
               multi.durations{event_index} = zeros(size(multi.onsets{event_index}));
               event_index = event_index +1;
            end   
            
            
       case 35
                       
            event_index = 1;
            
        
            
            if any(feedback_onsets)
                %all trials
                multi.names{event_index} = 'all_trials';
                multi.onsets{event_index} = feedback_onsets(~timeouts_latent_guess);
                multi.durations{event_index} = zeros(size(multi.onsets{event_index})); %impulse regressor (as opposed to boxcar regressor)
                event_index = event_index +1;
            end
            
            
%             if any(wins)
%             %win trials
%                 multi.names{event_index} = 'wins';
%                 multi.onsets{event_index} = feedback_onsets(wins);
%                 multi.durations{event_index} = zeros(size(multi.onsets{event_index})); %impulse regressor (as opposed to boxcar regressor)
%                 event_index = event_index +1;
%             end
% %             
%             if any(adversarial)
%             %adversarial trials
%                 multi.names{event_index} = 'adversarial';
%                 multi.onsets{event_index} = feedback_onsets(adversarial);
%                 multi.durations{event_index} = zeros(size(multi.onsets{event_index})); %impulse regressor (as opposed to boxcar regressor)
%                 event_index = event_index +1;
%             end
%             
%             if any(wins_adversarial)
%                 %win adversarial trials
%                 multi.names{event_index} = 'wins_adversarial';
%                 multi.onsets{event_index} = feedback_onsets(wins_adversarial);
%                 multi.durations{event_index} = zeros(size(multi.onsets{event_index})); %impulse regressor (as opposed to boxcar regressor)
%                 event_index = event_index +1;
%             end
%            

           multi.orth{1} = 0; % do not orthogonalise them            
                
            
            % nuisance @ reaction onset
            %
            multi.names{event_index} = 'reaction_onset';
            multi.onsets{event_index} = reaction_onsets;
            multi.durations{event_index} = zeros(size(multi.onsets{event_index}));
            event_index = event_index +1;
         

            % nuisance @ trial onset
            %
            multi.names{event_index} = 'trial_onset';
            multi.onsets{event_index} = trial_onsets;
            multi.durations{event_index} = zeros(size(multi.onsets{event_index}));
            event_index = event_index +1;
       

           % nuisance @ latent onset
            %
            multi.names{event_index} = 'latent_onset';
            multi.onsets{event_index} = latent_onsets;
            multi.durations{event_index} = zeros(size(multi.onsets{event_index}));
            event_index = event_index +1;
            
            % nuisance @ latent offset
            %
            multi.names{event_index} = 'latent_offset';
            multi.onsets{event_index} = latent_offsets;
            multi.durations{event_index} = zeros(size(multi.onsets{event_index}));
            event_index = event_index +1;
        
            
            if sum(timeouts_latent_guess) > 0
               multi.names{event_index} = 'feedback_onset_timeouts';
               multi.onsets{event_index} = feedback_onsets(timeouts_latent_guess); % timeouts only
               multi.durations{event_index} = zeros(size(multi.onsets{event_index}));
               event_index = event_index +1;
            end   
                               

            
             case 36
                       
            event_index = 1;
            
        
            
            if any(feedback_onsets)
                %all trials
                multi.names{event_index} = 'all_trials';
                multi.onsets{event_index} = feedback_onsets(~timeouts_latent_guess);
                multi.durations{event_index} = zeros(size(multi.onsets{event_index})); %impulse regressor (as opposed to boxcar regressor)
                event_index = event_index +1;
            end
            
            
%             if any(wins)
%             %win trials
%                 multi.names{event_index} = 'wins';
%                 multi.onsets{event_index} = feedback_onsets(wins);
%                 multi.durations{event_index} = zeros(size(multi.onsets{event_index})); %impulse regressor (as opposed to boxcar regressor)
%                 event_index = event_index +1;
%             end
% %             
%             if any(adversarial)
%             %adversarial trials
%                 multi.names{event_index} = 'adversarial';
%                 multi.onsets{event_index} = feedback_onsets(adversarial);
%                 multi.durations{event_index} = zeros(size(multi.onsets{event_index})); %impulse regressor (as opposed to boxcar regressor)
%                 event_index = event_index +1;
%             end
%             
            if any(wins_adversarial)
                %win adversarial trials
                multi.names{event_index} = 'wins_adversarial';
                multi.onsets{event_index} = feedback_onsets(wins_adversarial);
                multi.durations{event_index} = zeros(size(multi.onsets{event_index})); %impulse regressor (as opposed to boxcar regressor)
                event_index = event_index +1;
            end
           

           multi.orth{1} = 0; % do not orthogonalise them            
                
            
            % nuisance @ reaction onset
            %
            multi.names{event_index} = 'reaction_onset';
            multi.onsets{event_index} = reaction_onsets;
            multi.durations{event_index} = zeros(size(multi.onsets{event_index}));
            event_index = event_index +1;
         

            % nuisance @ trial onset
            %
            multi.names{event_index} = 'trial_onset';
            multi.onsets{event_index} = trial_onsets;
            multi.durations{event_index} = zeros(size(multi.onsets{event_index}));
            event_index = event_index +1;
       

           % nuisance @ latent onset
            %
            multi.names{event_index} = 'latent_onset';
            multi.onsets{event_index} = latent_onsets;
            multi.durations{event_index} = zeros(size(multi.onsets{event_index}));
            event_index = event_index +1;
            
            % nuisance @ latent offset
            %
            multi.names{event_index} = 'latent_offset';
            multi.onsets{event_index} = latent_offsets;
            multi.durations{event_index} = zeros(size(multi.onsets{event_index}));
            event_index = event_index +1;
        
            
            if sum(timeouts_latent_guess) > 0
               multi.names{event_index} = 'feedback_onset_timeouts';
               multi.onsets{event_index} = feedback_onsets(timeouts_latent_guess); % timeouts only
               multi.durations{event_index} = zeros(size(multi.onsets{event_index}));
               event_index = event_index +1;
            end   
              
            
            
           case 37
                       
            event_index = 1;
            
        
            
            if any(feedback_onsets)
                %all trials
                multi.names{event_index} = 'all_trials';
                multi.onsets{event_index} = feedback_onsets(~timeouts_latent_guess);
                multi.durations{event_index} = zeros(size(multi.onsets{event_index})); %impulse regressor (as opposed to boxcar regressor)
                event_index = event_index +1;
            end
            
            
            if any(wins)
            %win trials
                multi.names{event_index} = 'wins';
                multi.onsets{event_index} = feedback_onsets(wins);
                multi.durations{event_index} = zeros(size(multi.onsets{event_index})); %impulse regressor (as opposed to boxcar regressor)
                event_index = event_index +1;
            end
% %             
%             if any(adversarial)
%             %adversarial trials
%                 multi.names{event_index} = 'adversarial';
%                 multi.onsets{event_index} = feedback_onsets(adversarial);
%                 multi.durations{event_index} = zeros(size(multi.onsets{event_index})); %impulse regressor (as opposed to boxcar regressor)
%                 event_index = event_index +1;
%             end
%             
            if any(wins_adversarial)
                %win adversarial trials
                multi.names{event_index} = 'wins_adversarial';
                multi.onsets{event_index} = feedback_onsets(wins_adversarial);
                multi.durations{event_index} = zeros(size(multi.onsets{event_index})); %impulse regressor (as opposed to boxcar regressor)
                event_index = event_index +1;
            end
           

           multi.orth{1} = 0; % do not orthogonalise them            
                
            
            % nuisance @ reaction onset
            %
            multi.names{event_index} = 'reaction_onset';
            multi.onsets{event_index} = reaction_onsets;
            multi.durations{event_index} = zeros(size(multi.onsets{event_index}));
            event_index = event_index +1;
         

            % nuisance @ trial onset
            %
            multi.names{event_index} = 'trial_onset';
            multi.onsets{event_index} = trial_onsets;
            multi.durations{event_index} = zeros(size(multi.onsets{event_index}));
            event_index = event_index +1;
       

           % nuisance @ latent onset
            %
            multi.names{event_index} = 'latent_onset';
            multi.onsets{event_index} = latent_onsets;
            multi.durations{event_index} = zeros(size(multi.onsets{event_index}));
            event_index = event_index +1;
            
            % nuisance @ latent offset
            %
            multi.names{event_index} = 'latent_offset';
            multi.onsets{event_index} = latent_offsets;
            multi.durations{event_index} = zeros(size(multi.onsets{event_index}));
            event_index = event_index +1;
        
            
            if sum(timeouts_latent_guess) > 0
               multi.names{event_index} = 'feedback_onset_timeouts';
               multi.onsets{event_index} = feedback_onsets(timeouts_latent_guess); % timeouts only
               multi.durations{event_index} = zeros(size(multi.onsets{event_index}));
               event_index = event_index +1;
            end   
      
            
            
  % yes > no at feedback
           case 38 
            
            event_index = 1;
            
            save debug.mat; 
            
            assert(sum(yes) + sum(no) <= length(latent_onsets))
            assert(length(yes) == length(latent_onsets))
            assert(length(no) == length(latent_onsets))
            
             if any(yes)         
                % latent guess = yes
                multi.names{event_index} = ['yes'];
                multi.onsets{event_index} = feedback_onsets(yes);
                multi.durations{event_index} = zeros(size(multi.onsets{event_index})); %impulse regressor (as opposed to boxcar regressor)
                event_index = event_index +1;
             end
             
             if any(no)         
                % latent guess = no
                multi.names{event_index} = ['no'];
                multi.onsets{event_index} = feedback_onsets(no);
                multi.durations{event_index} = zeros(size(multi.onsets{event_index})); %impulse regressor (as opposed to boxcar regressor)
                event_index = event_index +1;
             end
             
             % nuisance @ reaction onset
            %
            multi.names{event_index} = 'reaction_onset';
            multi.onsets{event_index} = reaction_onsets;
            multi.durations{event_index} = zeros(size(multi.onsets{event_index}));
            event_index = event_index +1;

            % nuisance @ trial onset
            %
            multi.names{event_index} = 'trial_onset';
            multi.onsets{event_index} = trial_onsets;
            multi.durations{event_index} = zeros(size(multi.onsets{event_index}));
            event_index = event_index +1;
            
            % nuisance @ latent onset
            %
            multi.names{event_index} = 'latent_onset';
            multi.onsets{event_index} = latent_onsets;
            multi.durations{event_index} = zeros(size(multi.onsets{event_index}));
            
            % nuisance @ latent offset
            %
            multi.names{event_index} = 'latent_offset';
            multi.onsets{event_index} = latent_offsets;
            multi.durations{event_index} = zeros(size(multi.onsets{event_index}));
         
               
            if sum(timeouts_latent_guess) > 0
               multi.names{event_index} = 'feedback_onset_timeouts';
               multi.onsets{event_index} = feedback_onsets(timeouts_latent_guess); % timeouts only
               multi.durations{event_index} = zeros(size(multi.onsets{event_index}));
            end         
                         
            
         % RPE in adversarial and benevolent as continuous pmods --- figure
         % out!
        % 
        case 39 
            % 
            %
            event_index = 1;
             
                
            [RPE,psi] = get_latents(subj_original_idx, run);
            
            if any(adversarial & ~timeouts_latent_guess)
            %adversarial trials
                multi.names{event_index} = 'adversarial';
                multi.onsets{event_index} = feedback_onsets(adversarial & ~timeouts_latent_guess);
                multi.durations{event_index} = zeros(size(multi.onsets{event_index})); %impulse regressor (as opposed to boxcar regressor)
                multi.pmod(event_index).name{1} = 'RPE';
                multi.pmod(event_index).param{1} = RPE(adversarial & ~timeouts_latent_guess)';
                multi.pmod(event_index).poly{1} = 1; 
                event_index = event_index +1;
            end

            if any(benevolent & ~timeouts_latent_guess)
            %benevolent trials
                multi.names{event_index} = 'benevolent';
                multi.onsets{event_index} = feedback_onsets(benevolent & ~timeouts_latent_guess);
                multi.durations{event_index} = zeros(size(multi.onsets{event_index})); %impulse regressor (as opposed to boxcar regressor)
                multi.pmod(event_index).name{1} = 'RPE';
                multi.pmod(event_index).param{1} = RPE(benevolent & ~timeouts_latent_guess)';
                multi.pmod(event_index).poly{1} = 1; 
                event_index = event_index +1;
            end

           multi.orth{1} = 0; % do not orthogonalise them  
           %RPE = RPE(~timeouts_latent_guess
  

            % nuisance @ reaction onset
            %
            multi.names{event_index} = 'reaction_onset';
            multi.onsets{event_index} = reaction_onsets;
            multi.durations{event_index} = zeros(size(multi.onsets{event_index}));
            event_index = event_index +1;

            % nuisance @ feedback onset
            %
            multi.names{event_index} = 'trial_onset';
            multi.onsets{event_index} = trial_onsets;
            multi.durations{event_index} = zeros(size(multi.onsets{event_index}));
            event_index = event_index +1;
            
            % nuisance @ latent onset
            %
            multi.names{event_index} = 'latent_onset';
            multi.onsets{event_index} = latent_onsets;
            multi.durations{event_index} = zeros(size(multi.onsets{event_index}));
            event_index = event_index +1;
            
            % nuisance @ latent offset
            %
            multi.names{event_index} = 'latent_offset';
            multi.onsets{event_index} = latent_offsets;
            multi.durations{event_index} = zeros(size(multi.onsets{event_index}));
            event_index = event_index +1;
            
               
            if sum(timeouts_latent_guess) > 0
               multi.names{event_index} = 'feedback_onset_timeouts';
               multi.onsets{event_index} = feedback_onsets(timeouts_latent_guess); % timeouts only
               multi.durations{event_index} = zeros(size(multi.onsets{event_index}));
            end     
   
    % wins adversarial - wins benevolent
            case 40    
                       
            event_index = 1;
                      
            if any(wins & benevolent & ~timeouts_latent_guess)
                %win benevolent trials
                multi.names{event_index} = 'wins_benevolent';
                multi.onsets{event_index} = feedback_onsets(wins & benevolent & ~timeouts_latent_guess);
                multi.durations{event_index} = zeros(size(multi.onsets{event_index})); %impulse regressor (as opposed to boxcar regressor)
                event_index = event_index +1;
            end
            
            if any(wins & adversarial & ~timeouts_latent_guess)
                %win adversarial trials
                multi.names{event_index} = 'wins_adversarial';
                multi.onsets{event_index} = feedback_onsets(wins & adversarial & ~timeouts_latent_guess);
                multi.durations{event_index} = zeros(size(multi.onsets{event_index})); %impulse regressor (as opposed to boxcar regressor)
                event_index = event_index +1;
            end
            
            
                       multi.orth{1} = 0; % do not orthogonalise them            
                
            
            % nuisance @ reaction onset
            %
            multi.names{event_index} = 'reaction_onset';
            multi.onsets{event_index} = reaction_onsets;
            multi.durations{event_index} = zeros(size(multi.onsets{event_index}));
            event_index = event_index +1;
         

            % nuisance @ trial onset
            %
            multi.names{event_index} = 'trial_onset';
            multi.onsets{event_index} = trial_onsets;
            multi.durations{event_index} = zeros(size(multi.onsets{event_index}));
            event_index = event_index +1;
       

           % nuisance @ latent onset
            %
            multi.names{event_index} = 'latent_onset';
            multi.onsets{event_index} = latent_onsets;
            multi.durations{event_index} = zeros(size(multi.onsets{event_index}));
            event_index = event_index +1;
            
            % nuisance @ latent offset
            %
            multi.names{event_index} = 'latent_offset';
            multi.onsets{event_index} = latent_offsets;
            multi.durations{event_index} = zeros(size(multi.onsets{event_index}));
            event_index = event_index +1;
        
            
            if sum(timeouts_latent_guess) > 0
               multi.names{event_index} = 'feedback_onset_timeouts';
               multi.onsets{event_index} = feedback_onsets(timeouts_latent_guess); % timeouts only
               multi.durations{event_index} = zeros(size(multi.onsets{event_index}));
               event_index = event_index +1;
            end   
           
 % losses adversarial - losses benevolent
            case 41    
                       
            event_index = 1;
                      
            if any(losses & benevolent & ~timeouts_latent_guess)
                %win benevolent trials
                multi.names{event_index} = 'losses_benevolent';
                multi.onsets{event_index} = feedback_onsets(losses & benevolent & ~timeouts_latent_guess);
                multi.durations{event_index} = zeros(size(multi.onsets{event_index})); %impulse regressor (as opposed to boxcar regressor)
                event_index = event_index +1;
            end
            
            if any(losses & adversarial & ~timeouts_latent_guess)
                %win adversarial trials
                multi.names{event_index} = 'losses_adversarial';
                multi.onsets{event_index} = feedback_onsets(losses & adversarial & ~timeouts_latent_guess);
                multi.durations{event_index} = zeros(size(multi.onsets{event_index})); %impulse regressor (as opposed to boxcar regressor)
                event_index = event_index +1;
            end
           
           multi.orth{1} = 0; % do not orthogonalise them            
                
            
            % nuisance @ reaction onset
            %
            multi.names{event_index} = 'reaction_onset';
            multi.onsets{event_index} = reaction_onsets;
            multi.durations{event_index} = zeros(size(multi.onsets{event_index}));
            event_index = event_index +1;
         

            % nuisance @ trial onset
            %
            multi.names{event_index} = 'trial_onset';
            multi.onsets{event_index} = trial_onsets;
            multi.durations{event_index} = zeros(size(multi.onsets{event_index}));
            event_index = event_index +1;
       

           % nuisance @ latent onset
            %
            multi.names{event_index} = 'latent_onset';
            multi.onsets{event_index} = latent_onsets;
            multi.durations{event_index} = zeros(size(multi.onsets{event_index}));
            event_index = event_index +1;
            
            % nuisance @ latent offset
            %
            multi.names{event_index} = 'latent_offset';
            multi.onsets{event_index} = latent_offsets;
            multi.durations{event_index} = zeros(size(multi.onsets{event_index}));
            event_index = event_index +1;
        
            
            if sum(timeouts_latent_guess) > 0
               multi.names{event_index} = 'feedback_onset_timeouts';
               multi.onsets{event_index} = feedback_onsets(timeouts_latent_guess); % timeouts only
               multi.durations{event_index} = zeros(size(multi.onsets{event_index}));
               event_index = event_index +1;
            end   
            
            
            
     % wins adversarial - losses benevolent
            case 42    
                       
            event_index = 1;
                      
            if any(losses & benevolent & ~timeouts_latent_guess)
                %loss benevolent trials
                multi.names{event_index} = 'losses_benevolent';
                multi.onsets{event_index} = feedback_onsets(losses & benevolent & ~timeouts_latent_guess);
                multi.durations{event_index} = zeros(size(multi.onsets{event_index})); %impulse regressor (as opposed to boxcar regressor)
                event_index = event_index +1;
            end
            
            if any(wins & adversarial & ~timeouts_latent_guess)
                %win adversarial trials
                multi.names{event_index} = 'wins_adversarial';
                multi.onsets{event_index} = feedback_onsets(wins & adversarial & ~timeouts_latent_guess);
                multi.durations{event_index} = zeros(size(multi.onsets{event_index})); %impulse regressor (as opposed to boxcar regressor)
                event_index = event_index +1;
            end
             
                   multi.orth{1} = 0; % do not orthogonalise them            
                
            
            % nuisance @ reaction onset
            %
            multi.names{event_index} = 'reaction_onset';
            multi.onsets{event_index} = reaction_onsets;
            multi.durations{event_index} = zeros(size(multi.onsets{event_index}));
            event_index = event_index +1;
         

            % nuisance @ trial onset
            %
            multi.names{event_index} = 'trial_onset';
            multi.onsets{event_index} = trial_onsets;
            multi.durations{event_index} = zeros(size(multi.onsets{event_index}));
            event_index = event_index +1;
       

           % nuisance @ latent onset
            %
            multi.names{event_index} = 'latent_onset';
            multi.onsets{event_index} = latent_onsets;
            multi.durations{event_index} = zeros(size(multi.onsets{event_index}));
            event_index = event_index +1;
            
            % nuisance @ latent offset
            %
            multi.names{event_index} = 'latent_offset';
            multi.onsets{event_index} = latent_offsets;
            multi.durations{event_index} = zeros(size(multi.onsets{event_index}));
            event_index = event_index +1;
        
            
            if sum(timeouts_latent_guess) > 0
               multi.names{event_index} = 'feedback_onset_timeouts';
               multi.onsets{event_index} = feedback_onsets(timeouts_latent_guess); % timeouts only
               multi.durations{event_index} = zeros(size(multi.onsets{event_index}));
               event_index = event_index +1;
            end       
            
 
            
            
            
        % wins adversarial - wins benevolent w/out reaction onset regressor
            case 43    
                       
            event_index = 1;
                      
            if any(wins & benevolent & ~timeouts_latent_guess)
                %win benevolent trials
                multi.names{event_index} = 'wins_benevolent';
                multi.onsets{event_index} = feedback_onsets(wins & benevolent & ~timeouts_latent_guess);
                multi.durations{event_index} = zeros(size(multi.onsets{event_index})); %impulse regressor (as opposed to boxcar regressor)
                event_index = event_index +1;
            end
            
            if any(wins & adversarial & ~timeouts_latent_guess)
                %win adversarial trials
                multi.names{event_index} = 'wins_adversarial';
                multi.onsets{event_index} = feedback_onsets(wins & adversarial & ~timeouts_latent_guess);
                multi.durations{event_index} = zeros(size(multi.onsets{event_index})); %impulse regressor (as opposed to boxcar regressor)
                event_index = event_index +1;
            end
            
            
                       multi.orth{1} = 0; % do not orthogonalise them            
                
            
            % nuisance @ reaction onset
            %
%             multi.names{event_index} = 'reaction_onset';
%             multi.onsets{event_index} = reaction_onsets;
%             multi.durations{event_index} = zeros(size(multi.onsets{event_index}));
%             event_index = event_index +1;
%          

            % nuisance @ trial onset
            %
            multi.names{event_index} = 'trial_onset';
            multi.onsets{event_index} = trial_onsets;
            multi.durations{event_index} = zeros(size(multi.onsets{event_index}));
            event_index = event_index +1;
       

           % nuisance @ latent onset
            %
            multi.names{event_index} = 'latent_onset';
            multi.onsets{event_index} = latent_onsets;
            multi.durations{event_index} = zeros(size(multi.onsets{event_index}));
            event_index = event_index +1;
            
            % nuisance @ latent offset
            %
            multi.names{event_index} = 'latent_offset';
            multi.onsets{event_index} = latent_offsets;
            multi.durations{event_index} = zeros(size(multi.onsets{event_index}));
            event_index = event_index +1;
        
            
            if sum(timeouts_latent_guess) > 0
               multi.names{event_index} = 'feedback_onset_timeouts';
               multi.onsets{event_index} = feedback_onsets(timeouts_latent_guess); % timeouts only
               multi.durations{event_index} = zeros(size(multi.onsets{event_index}));
               event_index = event_index +1;
            end   
            
            
            
            
        otherwise
            assert(false, 'invalid glmodel -- should be one of the above');

    end % end of switch statement

   if save_output
       save('exploration_create_multi.mat'); % <-- DON'T DO IT! breaks on NCF... b/c of file permissions
   end
end

function [RPE,psi] = get_latents(subj_original_idx, run)
           load('model_output.mat', 'results', 'data');
           RPE = [];
           psi = [];
           for s = 1:length(data)
               if data(s).sub == subj_original_idx
                   w = data(s).run_num == run;
                   RPE = results(2).latents(s).rpe(w);
                   psi = results(2).latents(s).psi(w);
               end
           end
           assert(~isempty(RPE));
           save debug.mat;
           %assert(length(RPE) == length(multi.onsets{1}));
           %assert(length(psi) == length(multi.onsets{1}));
end










           
      % Win vs. Loss in Benevolent vs. Adversarial w/ psi PMOD
        % Condition @ feedback onset 
        % nuisance @ reaction and trial onset 
        % Contrast to run: 1) wins_benevolent - losses_benevolent, 2)
        % wins_adversarial - losses_adversarial, 3)wins_benevolent +
        % wins_adversarial - losses_benevolent - losses_adversarial (so all
        % wins - all losses), 4)wins_adversarial - losses_adversarial -
        % wins_benevolent + losses_benevolent
        
%         case 10 
%             
%             event_index = 1;
%                
%            
%             
%              if any(wins)         
%                 % wins for condition
%                 multi.names{event_index} = ['wins_',condition];
%                 multi.onsets{event_index} = feedback_onsets(wins);
%                 multi.durations{event_index} = zeros(size(multi.onsets{event_index})); %impulse regressor (as opposed to boxcar regressor)
%                 [RPE,psi] = get_latents(subj_original_idx, run);
%                 multi.pmod(event_index).name{1} = 'psi';
%                 multi.pmod(event_index).param{1} = psi(wins(~timeouts_latent_guess))';
%                 multi.pmod(event_index).poly{1} = 1; 
%                 event_index = event_index +1;
%              end
%            
%             if any(losses)
%                 % losses for condition
%                 multi.names{event_index} = ['losses_',condition];
%                 multi.onsets{event_index} = feedback_onsets(losses);
%                 multi.durations{event_index} = zeros(size(multi.onsets{event_index})); %impulse regressor (as opposed to boxcar regressor)
%                 [RPE,psi] = get_latents(subj_original_idx, run);
%                 multi.pmod(event_index).name{1} = 'psi';
%                 multi.pmod(event_index).param{1} = psi(losses(~timeouts_latent_guess))';
%                 multi.pmod(event_index).poly{1} = 1; 
%                 event_index = event_index +1;
%             end
% 
%             % nuisance @ reaction onset
%             %
%             multi.names{event_index} = 'reaction_onset';
%             multi.onsets{event_index} = reaction_onsets;
%             multi.durations{event_index} = zeros(size(multi.onsets{event_index}));
%             event_index = event_index +1;
% 
%             % nuisance @ trial onset
%             %
%             multi.names{event_index} = 'trial_onset';
%             multi.onsets{event_index} = trial_onsets;
%             multi.durations{event_index} = zeros(size(multi.onsets{event_index}));
%             event_index = event_index +1;
%          
%                
%             if sum(timeouts) > 0
%                multi.names{event_index} = 'feedback_onset_timeouts';
%                multi.onsets{event_index} = feedback_onsets(timeouts); % timeouts only
%                multi.durations{event_index} = zeros(size(multi.onsets{event_index}));
%             end
%              
%            
            % Win vs. Loss in Benevolent vs. Adversarial w/ psi PMOD
        % Condition @ feedback onset 
        % nuisance @ reaction and trial onset 
        % Contrast to run: 1) wins_benevolent - losses_benevolent, 2)
        % wins_adversarial - losses_adversarial, 3)wins_benevolent +
        % wins_adversarial - losses_benevolent - losses_adversarial (so all
        % wins - all losses), 4)wins_adversarial - losses_adversarial -
        % wins_benevolent + losses_benevolent
        
             
