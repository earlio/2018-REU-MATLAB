function lineage_matrix = build_lineage_matrix(sample_vector, total_time, generational_demographics, life_table, ...
    final_population, full_history_table)

    disp(generational_demographics)
    
    lineage_matrix = initialize_lineage_matrix(sample_vector, total_time);

    for current_time = (total_time):-1:2
        
        previous_generation = get_previous_population( generational_demographics, current_time);
        
        
        offspring_assignments = generate_parent_chances(life_table, previous_generation, generational_demographics, current_time)
        
        previous_generation{1} 
        
    end

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

end


% returns an array of cells of 2 by x matrices - EZ
function previous_generation = get_previous_population(generational_demographics, current_time)


    
    previous_demographics = generational_demographics( (current_time - 1), : )
    
    number_of_age_cohorts = length(previous_demographics)
    
    previous_generation = cell(number_of_age_cohorts);
    
    for cohort_age = 1:number_of_age_cohorts
        
        individual_id_list = 1:previous_demographics(cohort_age);
        
        individual_age_list = ( cohort_age - 1) * ones(1, previous_demographics(cohort_age));
        
        previous_generation{cohort_age} = [individual_id_list ; individual_age_list];
   
    end
    
end

function offspring_assignments = generate_parent_chances(life_table, previous_generation, generational_demographics, current_time)

    fecundity_by_age = life_table(:,3)
    
    how_many_currently_age_zero = generational_demographics(current_time, 1)
    offspring_assignments = zeros(2, how_many_currently_age_zero)

    number_of_age_cohorts = length(previous_generation);
    
    offspring_assignments_distributions = cell(number_of_age_cohorts);
    
    offspring_assignments_distributions_sum = 0;
    
    %  We randomly say something like 
    % "lineage #38 (age 5) has 6 offspring, #39 (age 5) has 2, etc." and 
    % check that the total # of offspring is large enough.
    while (offspring_assignments_distributions_sum < how_many_currently_age_zero)
    %TODO we probably want a loop count break in here in case fecundity is
    %so low that we never generate a sufficient number of offspring
        offspring_assignments_distributions_sum = 0;
        
        for i = 1:number_of_age_cohorts
        
            fprintf("age cohort : %f\n", i)
            offsprings_assignments_distributions{i} = poissrnd(fecundity_by_age(i), 1, length(previous_generation{i}(1,:)));
            
            offspring_assignments_distributions_sum = ...
                offspring_assignments_distributions_sum + sum(offsprings_assignments_distributions{i})
            
            celldisp(offsprings_assignments_distributions)
            %ceil(length(previous_generation{i}(1,:)) * fecundity_by_age(i))
        
        end
    end
    
end