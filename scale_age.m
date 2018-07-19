function [life_table]  = scale_age(life_table, number_of_years_to_shift) 

%   A function that pushes onset of adulthood (i.e. fecundity > 0) back
%   an arbitrary numbmer of years, while maintaining the overall shape of
%   the fecundity curve and boosting fecundities to keep the population 
%   growth rate unchanged. (For an explanation of our approach to scaling 
%   population growth, See Waples, Do, and Chopelet 2011; section 
%   'Development of the Approach', subsection 'Assumptions and notation')
    
%   OVERVIEW:
%   We transform a plot of fecundity (y-axis) and age (x-axis) into a
%   (piecewise linear) function (using the ‘fit’ function from the Curve 
%   Fitting Toolbox). Then we essentially squeeze that function to 
%   interpolate exactly the number of fecundities we want. (We could also 
%   us this technique to stretch the function so as to move onset-of-fecundity 
%   forward instead of backward.) Finally, we scale the life table growth 
%   rate back up to whatever it started out (for consistency's sake).

    original_fecundity_vector = life_table(:,3); 
    survival_vector = life_table(:,2);
    age_vector = life_table(:,1);
    count_of_age_cohorts = length(age_vector);
    
    survivorship_curve_vector = calc_survivorship_curve_vector(count_of_age_cohorts, survival_vector)
    
    raw_growth_rate = sum(original_fecundity_vector .* survivorship_curve_vector);
    
    fecundity_func = fit(age_vector,original_fecundity_vector,'linearinterp') ;

    % safe to delete
    % plot(f, age_vector, original_fecundity_vector);
    
    startValue = 1; % age 1
    endValue = count_of_age_cohorts;
    
     % shorter than the old length of adulthood
    new_length_of_adulthood = count_of_age_cohorts - number_of_years_to_shift;
    
    % new sampling points for interpolating new sample fecundties
    new_adult_fecundity_vector_sample_points = linspace(startValue, endValue, new_length_of_adulthood);
    
    % interpolated values for fecundities of the new adult age classes
    new_adult_fecundity_vector_values = fecundity_func(new_adult_fecundity_vector_sample_points);
    
    new_fecundity_vector = [zeros(number_of_years_to_shift,1); new_adult_fecundity_vector_values];
    
    recalculated_growth_rate = sum(new_fecundity_vector .* survivorship_curve_vector)
    
    % we scale up the new fecundity vector so that the new growth rate is 
    % the same as the old growth rate
    new_fecundity_vector = new_fecundity_vector .* (raw_growth_rate/recalculated_growth_rate);

   % safe to delete; checks that new growth rate is the same as old
   % final_growth_rate = sum(new_fecundity_vector .* survivorship_curve_vector);

   life_table(:,3) = new_fecundity_vector;

end

% This formula produces the l-sub-x vector descrbied in Waples 2011)
function [survivorship_curve_vector] = calc_survivorship_curve_vector(count_of_age_cohorts, survival_vector);

    survivorship_curve_vector = ones(count_of_age_cohorts, 1);
    
    for i = 2:count_of_age_cohorts
        survivorship_curve_vector(i) = survivorship_curve_vector(i - 1) * survival_vector(i - 1);
    end
end