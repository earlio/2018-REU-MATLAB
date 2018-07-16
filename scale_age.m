function [new_life_table]  = scale_age(original_life_table, number_of_years_to_shift) 

%   First we calcule population growth rate so we can restore it
%   (See Waples, Do, and Chopelet 2011; section 'Development of the
%   Approach', subsection 'Assumptions and notation')
    
    original_fecundity_vector = original_life_table(:,3); 
    original_survival_vector = original_life_table(:,2);
 
    original_age_vector = original_life_table(:,1);
    
    count_of_age_cohorts = length(original_age_vector);
    
    survivorship_curve_vector = calc_survivorship_curve_vector(count_of_age_cohorts, original_survival_vector)
    
%   First we calcule population growth rate so we can restore it
%   (See Waples, Do, and Chopelet 2011; section 'Development of the
%   Approach', subsection 'Assumptions and notation')
    
    raw_growth_rate = sum(original_fecundity_vector .* survivorship_curve_vector);
    
    f = fit(original_age_vector,original_fecundity_vector,'linearinterp') ;
    

    
    
    plot(original_age_vector, original_fecundity_vector);
    

    

    
    startValue = 1;
    endValue = count_of_age_cohorts;
    new_length_of_adulthood = count_of_age_cohorts - number_of_years_to_shift; % shorter than old length ofadulthood
    new_adult_fecundity_vector_sample_points = linspace(startValue, endValue, new_length_of_adulthood);
    new_adult_fecundity_vector_values = f(new_adult_fecundity_vector_sample_points);
    new_fecundity_vector = [zeros(number_of_years_to_shift,1); new_adult_fecundity_vector_values];
    
    recalculated_growth_rate = sum(new_fecundity_vector .* survivorship_curve_vector)
    
    
    new_life_table = original_life_table;
    
   % old_fecundity_vector = original_life_table

% extract a formula from the fecundity vector
%%% f=fit(cdate,pop,'poly2') ; plot(f,cdate,pop)

% use formula to populate new_fecundity_vector

% multiply each element of new_fecundity_vector by 
% sum_old_vec / sum_new_vec

% install new_fecundity_vector

end

function [survivorship_curve_vector] = calc_survivorship_curve_vector(count_of_age_cohorts, survival_vector);

    survivorship_curve_vector = ones(count_of_age_cohorts, 1);
    
    for i = 2:count_of_age_cohorts
        survivorship_curve_vector(i) = survivorship_curve_vector(i - 1) * survival_vector(i - 1);
    end
end