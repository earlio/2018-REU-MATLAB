% scales life table growth using the method described in Waples 2011
function [original_scaling_factor, life_table] = life_table_scale(life_table, scaling_factor)

    original_fecundity_vector = life_table(:,3); 
    survival_vector = life_table(:,2);

    count_of_age_cohorts = length(survival_vector);
    
    survivorship_curve_vector = calc_survivorship_curve_vector(count_of_age_cohorts, survival_vector);
    
    original_scaling_factor = sum(original_fecundity_vector .* survivorship_curve_vector);
    
    new_fecundity_vector = original_fecundity_vector .* (scaling_factor/original_scaling_factor);
        
    life_table(:,3) = new_fecundity_vector;
end

% Produces the L-sub-x vector descrbied in Waples 2011
function [survivorship_curve_vector] = calc_survivorship_curve_vector(count_of_age_cohorts, survival_vector);

    survivorship_curve_vector = ones(count_of_age_cohorts, 1);
    
    for i = 2:count_of_age_cohorts
        survivorship_curve_vector(i) = survivorship_curve_vector(i - 1) * survival_vector(i - 1);
    end
end