function lineage_matrix = build_lineage_matrix(k_vector, current_population, ...
    life_table, total_time, generational_demographics, full_history_table)
%We return a t-by-k matrix of k samples, also known as lineages. We
%number each lineage and as we go backward in time, more and more lineages
%will share the same number.
    
    lineage_matrix = initialize_lineage_matrix(k_vector, total_time);
    
    % EZ number of columns in our lineage matrix, also known as k.
    k_size = size(lineage_matrix, 2);
    
    for current_time = (total_time):-1:2
        
        disp(lineage_matrix(current_time,:,1))
        disp("~~~")
        disp(lineage_matrix(current_time,:,2))
        
        for i = 1:k_size
        % first we age everyone backward one step and copy thier lineage
     
            lineage_matrix(current_time - 1, i, 2) = (lineage_matrix(current_time, i, 2) - 1);
            % we assume initially that there won't be a coalescent event
            lineage_matrix(current_time - 1, i, 1) = lineage_matrix(current_time, i, 1);
        end
        
        parent_chances_by_age = generate_parent_chances(full_history_table, ...
            life_table, current_time)
        
        for i = 1:k_size
            
            if (lineage_matrix(current_time - 1, i, 2) == -1)
                
                disp(lineage_matrix(:,:,1));
                disp(lineage_matrix(:,:,2));
                
                winning_ticket = randi(size(parent_chances_by_age, 2));

                sampled_age = parent_chances_by_age(winning_ticket)
                % remove winning ticket from pool of tickets
                parent_chances_by_age(winning_ticket) = [];
                

               
                %% determine whether parent is in sample or outside of sample %%
                
                % get count of all population at time t - 1
                previous_year_population = full_history_table(current_time - 1, :)
                population_age_cohort_size = size(previous_year_population(previous_year_population == sampled_age), 2)
                
                previous_year_sample_ages = lineage_matrix(current_time - 1, :,2)
                sample_age_cohort_size = size(previous_year_sample_ages(previous_year_sample_ages == sampled_age), 2)
                
                raffle_ticket_parent_is_in_sample = randi(population_age_cohort_size)
                
                is_in_sample = false;
                if (raffle_ticket_parent_is_in_sample  <= sample_age_cohort_size)
                    is_in_sample = true
                else
                    is_in_sample = false
                end
  
                if (is_in_sample)
   
                    potential_parents = lineage_matrix(current_time -1, :, :)
                    potential_parent_lineages = potential_parents(potential_parents(:,:,2) == sampled_age)
                    
                    count_of_potential_parents = size(potential_parent_lineages, 2)
                    
                    random_number_choose_parent = randi(count_of_potential_parents) 
                    
                    parent = potential_parent_lineages(random_number_choose_parent)
                    
                    lineage_matrix(current_time -1, i, 1) = parent
   
                end
                
                % We musn't forget to update the lineage with it's new age
                lineage_matrix(current_time - 1, i, 2) = sampled_age;

            end
           
            
     
        end
        
        % now we drop raffle tickets into a hat to determine which lineage
        % currently of reproductive age (wether in the sample of k lineages
        % or outside the sample) gets to be the ancestor of the sample
        % lineages that have aged back to -1
        
        % then, for each -1, we select a parent age and select an ancestor.
        % If the ancestor is alive in one of our sample lineages, then a
        % COALESCENT EVENT has occured.
        
        % check for and log coalescent events.
        end
    
    disp(lineage_matrix(:,:,1))
    disp(lineage_matrix(:,:,2))
end

% EZ for now we'll just look up the previous population rather than
% calculate it

function prev_population = build_previous_population(full_hist_table, tot_time, curr_time)
    prev_population  = full_hist_table(tot_time - curr_time, :);
end
    
function age_of_parent = calculate_age_of_parent(previous_population, life_table)

    age_of_parent = 1;

    % generate an array with total number of offspring for each age cohort,
    % but summed, so that if age_1 = 20 and age_2 = 30 then age_1 = 20 and
    % age_2 = 50
    
    % generate a random number between 1 and the maximum value in the
    % array above
    
    % loop through the above array with a counter and test wether our
    % random mumber is greater than the summed number of offspring
end
    
function lineage_matrix = initialize_lineage_matrix(k_vector, total_time)
% for now we're hard-coding the initial vector

    lineage_matrix = NaN(total_time, 3);
    lineage_matrix(:,:,2) = NaN;

        lineage_matrix(end, 1, 1) = 1;
        lineage_matrix(end, 1, 2) = 0;
        lineage_matrix(end, 2, 1) = 2;
        lineage_matrix(end, 2, 2) = 0;        
        lineage_matrix(end, 3, 1) = 3;
        lineage_matrix(end, 3, 2) = 2;
        lineage_matrix(:,:,1);
        lineage_matrix(:,:,2);
    
%     age = 0; % age of lineage to find
% 
%     found_one = 0;
%     while found_one == 0
% 
%         index = randi(total_pop_sizes(numgen,1)); % get random integer 1:total population size at time numgen
%         if pop_in_present(1,index) == age % lineage must be given age 
%             ancestors(numgen,1) = index;
%             found_one = 1;
%         else
%             found_one = 0; % keep looking for age zero individual
%         end
%     end % while
% 
%     % sample a second and different lineage in the present
%     age = 0; % age of lineage to find
% 
%     found_one = 0;
%     while found_one == 0
% 
%         index = randi(total_pop_sizes(numgen,1)); % get random integer 1:total population size at time numgen
%         if (pop_in_present(1,index) == age) && (index ~= ancestors(numgen,1)) % lineage must be given age and also not same as already sampled
%             ancestors(numgen,2) = index;
%             found_one = 1;
%         else
%             found_one = 0; % keep looking
%         end
%     end % while

end

function parent_chances = generate_parent_chances(full_history_table, life_table, current_time)
        
        parent_chances = [];
        previous_year_population = full_history_table(current_time - 1, :);
        
        % for each age cohort ...
        for i = 1:size(life_table, 1)
             
            % count the number of lineages ...
            cohort = size(previous_year_population(previous_year_population == i - 1), 2);
            
            % and multiply by fecundity ... 
            potential_offspring_count_by_cohort = ceil(cohort * life_table(i, 3));
            
            % then add one 'ticket', with the age of the current age
            % current age chort, for each potential_offspring
            for j = 1:potential_offspring_count_by_cohort
                
                parent_chances = [parent_chances, i - 1];
                
                
                
            end
            
            
            
        end
        
        
        
end