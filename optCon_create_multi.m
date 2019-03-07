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
    %ITI onset ??? ask if we should add
    wins = strcmp(data.feedback(which_trials),'1')'; %win trials
    losses = strcmp(data.feedback(which_trials), '0')'; %loss trials
    timeouts = strcmp(data.feedback(which_trials), 'NA')'; %missing trials
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

            % nuisance @ feedback onset
            %
            multi.names{3} = 'trial_onset';
            multi.onsets{3} = trial_onsets;
            multi.durations{3} = zeros(size(multi.onsets{3}));
           
            
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
           
           [RPE,psi] = get_latents(subj_original_idx, run, multi);
           
           multi.pmod(1).name{1} = 'RPEpsi';
           multi.pmod(1).param{1} = RPE' .* psi';
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
            
               
            if sum(timeouts_latent_guess) > 0
               multi.names{4} = 'feedback_onset_timeouts';
               multi.onsets{4} = feedback_onsets(timeouts_latent_guess); % timeouts only
               multi.durations{4} = zeros(size(multi.onsets{4}));
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
           
           [RPE,psi] = get_latents(subj_original_idx, run, multi);
                  
          
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
            
               
            if sum(timeouts_latent_guess) > 0
               multi.names{4} = 'feedback_onset_timeouts';
               multi.onsets{4} = feedback_onsets(timeouts_latent_guess); % timeouts only
               multi.durations{4} = zeros(size(multi.onsets{4}));
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
           
           [RPE,psi] = get_latents(subj_original_idx, run, multi);

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
            
               
            if sum(timeouts_latent_guess) > 0
               multi.names{4} = 'feedback_onset_timeouts';
               multi.onsets{4} = feedback_onsets(timeouts_latent_guess); % timeouts only
               multi.durations{4} = zeros(size(multi.onsets{4}));
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
           
           [RPE,psi] = get_latents(subj_original_idx, run, multi);
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
           disp(psi);
           disp(VS);
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
            
               
            if sum(timeouts_latent_guess) > 0
               multi.names{4} = 'feedback_onset_timeouts';
               multi.onsets{4} = feedback_onsets(timeouts_latent_guess); % timeouts only
               multi.durations{4} = zeros(size(multi.onsets{4}));
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
           
           [RPE,psi] = get_latents(subj_original_idx, run, multi);
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
            
               
            if sum(timeouts_latent_guess) > 0
               multi.names{4} = 'feedback_onset_timeouts';
               multi.onsets{4} = feedback_onsets(timeouts_latent_guess); % timeouts only
               multi.durations{4} = zeros(size(multi.onsets{4}));
            end         
            
            
           
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
               

        otherwise
            assert(false, 'invalid glmodel -- should be one of the above');

    end % end of switch statement

   if save_output
       save('exploration_create_multi.mat'); % <-- DON'T DO IT! breaks on NCF... b/c of file permissions
   end
end

function [RPE,psi] = get_latents(subj_original_idx, run, multi)
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
           assert(length(RPE) == length(multi.onsets{1}));
           assert(length(psi) == length(multi.onsets{1}));
end

