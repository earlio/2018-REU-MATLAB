function time_to_mrca = calculate_mrca(age_a, age_b, current_population, life_table, total_time, generational_demographics, full_history_table)

    time_to_mrca = -1;
    
    
    for t = 1:total_time
        
        % EZ we may want previous population to persist throughout this loop
        previous_population = []
            
        % EZ if age_a = 0 we know we've got to assign a parent
        if (age_a == 0)
            
            % EZ we need a a probability table for offspiring
            
            previous_population = build_previous_population(full_history_table, total_time, t)
            
            age_of_parent = calculate_age_of_parent(previous_population, life_table)
            
        end
    end
end

% EZ for now we'll just look up the previous population rather than
% calculate it

function previous_population = build_previous_population(full_history_table, total_time, current_time)
    previous_population  = full_history_table(total_time - current_time, :)
end
    
function age_of_parent = calculate_age_of_parent(previous_population, life_table)

    % generate an array with total number of offspring for each age cohort,
    % but summed, so that if age_1 = 20 and age_2 = 30 then age_1 = 20 and
    % age_2 = 50
    
    % generate a random number between 1 and the maximum value in the
    % array above
    
    % loop through the above array with a counter and test wether our
    % random mumber is greater than the summed number of offspring
end
    