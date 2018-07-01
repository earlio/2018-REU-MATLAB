function lineage_matrix = build_lineage_matrix(sample_vector, total_time, generational_demographics, life_table, ...
    final_population, full_history_table)

    disp(generational_demographics)
    
    lineage_matrix = initialize_lineage_matrix(sample_vector, total_time);

    for current_time = (total_time):-1:2
        
        previous_generation = get_previous_population( generational_demographics, current_time);
        
        
        parent_chances = generate_parent_chances(life_table, previous_generation, generational_demographics, current_time)
        
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

function parent_chances = generate_parent_chances(life_table, previous_generation, generational_demographics, current_time)

    fecundity_vector = life_table(:,3)
    
    number_of_age_cohorts = length(previous_generation)
    
    parent_chance_count_by_cohort = zeros(1, size(life_table, 1));
    
    parent_chances = zeros(2, generational_demographics(current_time -1, 1);
    
    % # of rows in life_table should be == length(previous_generation)
    for i = 1:number_of_age_cohorts
        
        parent_chances_by_cohort(i) = ceil(length(previous_generation{i}(1,:)) * fecundity_vector(i))
        
    end
end