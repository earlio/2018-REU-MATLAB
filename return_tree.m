function lineage_matrix = build_lineage_matrix(k_vector, current_population, ...
    life_table, total_time, generational_demographics, full_history_table)
%We're returning a t-by-k matrix of k samples, also known as lineages. We
%number each lineage and as we go backward in time, more and more lineages
%will share the same number.
    
    lineage_matrix = initialize_lineage_matrix(k_vector, total_time);
    
    % EZ number of columns in our lineage matrix, also known as k.
    k_size = size(lineage_matrix, 2);
    
    for current_time = (total_time):-1:2
        
        disp(lineage_matrix(current_time,:,1))
        disp(lineage_matrix(current_time,:,2))
        
        for i = 1:k_size
        % first we age everyone backward one step and copy thier lineage
     
            lineage_matrix(current_time - 1, i, 2) = (lineage_matrix(current_time, i, 2) - 1);
            lineage_matrix(current_time - 1, i, 1) = lineage_matrix(current_time, i, 1);
        end

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
    
% returns a total_time x 'sample size' x 2 matrix of -1's, with the sample
% demographics as the bottom row.
function lineage_matrix = initialize_lineage_matrix(k_vector, total_time)
        % for now we'll just return for k = 0, 0

         % EZ here is where we'd parse k_vector so we don't have to hard-code '2'
    % here
    lineage_matrix = NaN(total_time, 2);
    lineage_matrix(:,:,2) = NaN;

        lineage_matrix(end, 1, 1) = 1;
        lineage_matrix(end, 1, 2) = 0;
        lineage_matrix(end, 2, 1) = 2;
        lineage_matrix(end, 2, 2) = 0;
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