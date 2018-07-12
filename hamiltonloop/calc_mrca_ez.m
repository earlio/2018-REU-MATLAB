
% INPUTS:
% genealogy_m - contains number and ages of sampled lineages
% life table - contains fecundity data
% generational_demographics - contains population data
function  [mrca, complete_genealogy, coal_events] = calc_mrca_ez(genealogy_m, life_table, generational_demographics)

    


    % we'll make total_time = 5000 later
    
% START modifications required to work with calling program Main_mbh.m
    generational_demographics = transpose(generational_demographics);
    maximum_age = size(life_table, 1);
    % this function was designed for life tables starting at age 0 ...
    life_table(:, 1) = life_table(:, 1) - 1;
    
    % since we're starting with k = 2 I can manually build dample_vector
    sample_vector = zeros(1, maximum_age);
    first = genealogy_m(end,1,2) + 1;
    second = genealogy_m(end,2,2) + 1;
    sample_vector(first) = sample_vector(first) + 1;
    sample_vector(second) = sample_vector(second) + 1;
    
    total_time = size(genealogy_m,1);
    mrca = total_time;
    
    % the calling program isn't currently using these values
    complete_genealogy = 1;
    coal_events = 1;
    
% END modifications


    lineage_matrix = initialize_lineage_matrix(sample_vector, total_time);
    sample_size = length(lineage_matrix(1,:,1));

    for current_time = (total_time):-1:2
        
        previous_generation = interpolate_previous_population( generational_demographics, current_time);
        offspring_assignments = generate_parent_chances(life_table, previous_generation, lineage_matrix, current_time);
       
       
       for k = 1:sample_size
   
           age_of_sample = lineage_matrix(current_time, k, 2);
           
           if (age_of_sample > 0)
                
              %note: age_of_sample == 3 means we'd be looking up cohort 3+1
              %IF we wanted age 3, but we want age 2 so we look up cohort 3
              number_of_cohort_members_remaining = length(previous_generation{age_of_sample}(:));
              random_index = randi(number_of_cohort_members_remaining);
              
              lineage_matrix(current_time - 1, k, 1) = lineage_matrix(current_time, k, 1);
              lineage_matrix(current_time - 1, k, 2) = age_of_sample - 1;
              % reminder: age = z implies age group = z + 1
              lineage_matrix(current_time - 1, k, 3) = previous_generation{age_of_sample}(random_index); 
              
              % reminder: ages are from 0 to x but our array is 1 to x+1
              previous_generation{age_of_sample}(random_index) = [];
            
           end
           
           % NaN == NaN evaluates to false, so we need to use isnan()
           if (isnan(age_of_sample)) % do nothing; a merger has already taken place
             %  lineage_matrix(current_time - 1, k, 1) = lineage_matrix(current_time, k, 1);
               lineage_matrix(current_time - 1, k, 1) = NaN;
               lineage_matrix(current_time - 1, k, 2) = NaN;
               lineage_matrix(current_time - 1, k, 3) = NaN;
           end 
           
       end
       
       for k = 1:sample_size
           age_of_sample = lineage_matrix(current_time, k, 2);
           
           
           if (age_of_sample == 0)
               
               count_of_available_parent_assignments = length(offspring_assignments(1,:));
               random_parent_assignment_index = randi(count_of_available_parent_assignments);
               
               % actual parental age; age could be zero if fecundity allows
               parental_age = offspring_assignments(1,random_parent_assignment_index);
               parental_index = offspring_assignments(2,random_parent_assignment_index);
               offspring_assignments(:,random_parent_assignment_index) = [];
               
               % if the new parent is already in the sample, coalescent!
               
               is_coalescent = false;
               
               % O(k^2) :( maybe speed this up via a faster equality check?
               for potential_match = 1:sample_size
                   if ((lineage_matrix(current_time - 1, potential_match, 2) == parental_age) && ...
                       (lineage_matrix(current_time - 1, potential_match, 3) == parental_index ))
                      
                        % the new lineage of k is the lineage of the match
                        lineage_matrix( current_time - 1, k, 1 ) = lineage_matrix(current_time - 1, potential_match, 1);
                        lineage_matrix( current_time - 1, k, 2 ) = NaN;
                        lineage_matrix( current_time - 1, k, 3 ) = NaN;
                        
                        is_coalescent = true;
                        break
                        
                   end
               end
               
               
               % in the new parent is not
               
               if (is_coalescent == false)
                   
                   lineage_matrix(current_time - 1, k, 1) = lineage_matrix(current_time, k, 1);
                   lineage_matrix(current_time - 1, k, 2) = parental_age;
                   lineage_matrix(current_time - 1, k, 3) = parental_index;
               end
           end
       end
       
       % do we have an mrca?
       % rewrite this by using a counter after is_coalescent == true
       raw_lineage_numbers = (lineage_matrix(current_time -1, :, 1));
       actual_lineage_numbers = raw_lineage_numbers(isfinite(raw_lineage_numbers));

       if (length(actual_lineage_numbers) == 1)
           mrca = total_time - current_time + 1;
           break;
       end
       
      
    end
    
   % disp(lineage_matrix)
    
end

   
function lineage_matrix = initialize_lineage_matrix(k_vector, total_time)
% time x sample-size x 3; (sample_lineage, age, position in age cohort)

% for now we're hard-coding the initial vector

     sample_size = sum(k_vector);
     lineage_matrix = -1 * ones(total_time, sample_size ,3);

    lineage_id = 1;
    for i = 1:length(k_vector)
        for j = 1:k_vector(i)
            lineage_matrix(end, lineage_id, 1) = lineage_id;
            % we're setting age of sample here
            lineage_matrix(end, lineage_id, 2) = i - 1;
            lineage_id = lineage_id + 1;
            
        end
    end


end


% returns an array of cells of 2 by x matrices - EZ
function previous_generation = interpolate_previous_population(generational_demographics, current_time)
    
    previous_demographics = generational_demographics( (current_time - 1), : );
    
    number_of_age_cohorts = length(previous_demographics);
    
    previous_generation = cell(number_of_age_cohorts, 1);
    
    for cohort_age = 1:number_of_age_cohorts
        
        individual_id_list = 1:previous_demographics(cohort_age);
        
        % deprecated
        % individual_age_list = ( cohort_age - 1) * ones(1, previous_demographics(cohort_age));
        
        previous_generation{cohort_age} = individual_id_list;
   
    end
    
end

function offspring_assignments = generate_parent_chances(life_table, previous_generation, lineage_matrix, current_time)
    % returns 2 rows: row 1 = actual age, second row = 1st, 2nd, 4th of
    % that age; each column is a 'raffle ticket' in Darwin's lottery.
    fecundity_by_age = life_table(:,3);
    
   %  the number of age_zero in the sample, 
   % TODO non-critical not sure why debugger shows lineage_submatrix as zero
   lineage_submatrix = lineage_matrix(lineage_matrix(current_time, :, 2) == 0);
   count_of_currently_age_zero = length(lineage_submatrix(1,:,1));
 %   count_of_currently_age_zero_in_sample = 1

    number_of_age_cohorts = length(previous_generation);
    
    % "[age 5, lineage 2] has 6 offspring, [age 5, lineage 3] has 2, etc." 
    offspring_assignments_distributions = cell(number_of_age_cohorts, 1);
    
    offspring_assignments_distributions_sum = 0;

    % check that the total # of offspring is large enough.
    attemptCounter = 1;
    while (offspring_assignments_distributions_sum < count_of_currently_age_zero)
        
   %     fprintf("attempt #%d\n", attemptCounter)
        if (attemptCounter > 100)
            disp("too many attempts!")
            break
        end
        attemptCounter = attemptCounter + 1;
    %TODO we want a loop count break in here in case fecundity is very low
    
        offspring_assignments_distributions_sum = 0;
        
        for i = 1:number_of_age_cohorts
        
          %  fprintf("age cohort : %f\n", i)
            offspring_assignments_distributions{i} = poissrnd(fecundity_by_age(i), 1, length(previous_generation{i}(1,:)));
            
            offspring_assignments_distributions_sum = ...
                offspring_assignments_distributions_sum + sum(offspring_assignments_distributions{i});
            
           % celldisp(offspring_assignments_distributions);
            %ceil(length(previous_generation{i}(1,:)) * fecundity_by_age(i))
        
        end
    end
    
    offspring_assignments = zeros(2, offspring_assignments_distributions_sum);
    offspring_assignments_counter = 1;
    
    for age_cohort = 1:number_of_age_cohorts
        
        number_of_lineages = length(offspring_assignments_distributions{age_cohort}(:));
        for lineage = 1:number_of_lineages
            
            individual_number_of_offspring = offspring_assignments_distributions{age_cohort}(1,lineage);
            
            for offspring_index = 1:individual_number_of_offspring
                % remember age cohorts are 1 to z but ages are 0 to z - 1
                offspring_assignments(1:2, offspring_assignments_counter) = [age_cohort - 1, lineage];
                offspring_assignments_counter = offspring_assignments_counter + 1;
                
            end
        end
        
        
    end
    
    % don't think we need this
%     while (length(offspring_assignments(1,:)) > count_of_currently_age_zero)
%       %  disp("too many offspring ...")
%         remove_index = randi(length(offspring_assignments(1,:)));
%         offspring_assignments( :, remove_index) = [];
%     end
    
    
end