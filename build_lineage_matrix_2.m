function lineage_matrix = build_lineage_matrix(sample_vector, total_time, generational_demographics, life_table, ...
    final_population, full_history_table)

    disp(generational_demographics)
    
    lineage_matrix = initialize_lineage_matrix(sample_vector, total_time)
    sample_size = length(lineage_matrix(1,:,1))

    for current_time = (total_time):-1:2
        
        previous_generation = get_previous_population( generational_demographics, current_time);
        celldisp(previous_generation);
        
        disp("prev gen run")
        
        offspring_assignments = generate_parent_chances(life_table, previous_generation, generational_demographics, current_time)
        
        previous_generation{1} 
        
        
        
        % for each each k in the sample
        
            % if age > 0
                % pop an index at random from previous_generation
                
            % if age = 0
                %pop an index at random from offspring_assignments
        
            % if age < 0
                % implies a merger has already taken place: pass same age
                % and lineage backward
        
       %
       
       disp(sample_size)
       for k = 1:sample_size
           disp(k)
           disp("lineage matrix(current_time, k, 1)")
           lineage_matrix(current_time, k, 1)
           age_of_sample = lineage_matrix(current_time, k, 2)
           if (age_of_sample > 0)
 
              index = randi(previous_generation{age_of_sample + 1}(1,:))
             %  lineage_matrix(current_time
           end
           
           
       end
      % lineage_matrix
    end

end

   
function lineage_matrix = initialize_lineage_matrix(k_vector, total_time)
% time by sample-size by 3; (sample_lineage, age, nth-element by age group)

% for now we're hard-coding the initial vector

    sample_size = 3
    lineage_matrix = NaN(total_time, sample_size ,3);
   % lineage_matrix(:,:,2) = NaN;

        lineage_matrix(end, 1, 1) = 1;
        lineage_matrix(end, 1, 2) = 0;
        lineage_matrix(end, 2, 1) = 2;
        lineage_matrix(end, 2, 2) = 0;        
        lineage_matrix(end, 3, 1) = 3;
        lineage_matrix(end, 3, 2) = 2;
        lineage_matrix(:,:,1);
        lineage_matrix(:,:,2);

end


% returns an array of cells of 2 by x matrices - EZ
function previous_generation = get_previous_population(generational_demographics, current_time)


    
    previous_demographics = generational_demographics( (current_time - 1), : );
    
    number_of_age_cohorts = length(previous_demographics);
    
    previous_generation = cell(number_of_age_cohorts);
    
    for cohort_age = 1:number_of_age_cohorts
        
        individual_id_list = 1:previous_demographics(cohort_age);
        
        individual_age_list = ( cohort_age - 1) * ones(1, previous_demographics(cohort_age));
        
        previous_generation{cohort_age} = [individual_id_list ; individual_age_list];
   
    end
    
end

function offspring_assignments = generate_parent_chances(life_table, previous_generation, generational_demographics, current_time)

    fecundity_by_age = life_table(:,3);
    
    count_of_currently_age_zero = generational_demographics(current_time, 1);
    

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
        
        number_of_lineages = length(offspring_assignments_distributions{age_cohort}(1,:));
        for lineage = 1:number_of_lineages
            
            individual_number_of_offspring = offspring_assignments_distributions{age_cohort}(1,lineage);
            
            for offspring_index = 1:individual_number_of_offspring
                offspring_assignments(1:2, offspring_assignments_counter) = [age_cohort, lineage];
               % disp("just made an assignment... ")
                offspring_assignments_counter = offspring_assignments_counter + 1;
                
            end
        end
        
        
    end
    
    while (length(offspring_assignments(1,:)) > count_of_currently_age_zero)
      %  disp("too many offspring ...")
        remove_index = randi(length(offspring_assignments(1,:)));
        offspring_assignments( :, remove_index) = [];
    end
    
    
end