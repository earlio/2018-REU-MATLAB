% Standalone program that takes the name of a life table as input and
% writes several alternative life tables (.xlsx) to disk with onset of
% adulthood scaled back to various levels.

file_path_name = 'Sample_Pine_life_table';
extension = '.xlsx';
full_file_path_name = strcat(file_path_name, extension);

file_contents = importdata(full_file_path_name);
life_table = file_contents.data;

age_vector = life_table(:,1);
count_of_age_cohorts = length(age_vector);
shift_point5 = floor(count_of_age_cohorts * 0.5);
shift_point75 = floor(count_of_age_cohorts * 0.75);
shift_point825 = floor(count_of_age_cohorts * 0.825);
shift_point9 = floor(count_of_age_cohorts * 0.9);

if (shift_point5 == 0) 
    disp("error: not enough age cohorts to push onset of fecundity back by 50% of life span. Exiting program.");
    return;
end

file_name_to_print =  strcat(file_path_name, "_shift_point5", extension);
make_shifted_table(life_table, shift_point5, file_name_to_print);

file_name_to_print =  strcat(file_path_name, "_shift_point75", extension);
make_shifted_table(life_table, shift_point75, file_name_to_print);

file_name_to_print =  strcat(file_path_name, "_shift_point825", extension);
make_shifted_table(life_table, shift_point825, file_name_to_print);

file_name_to_print =  strcat(file_path_name, "_shift_point9", extension);
make_shifted_table(life_table, shift_point9, file_name_to_print);


function [] = make_shifted_table(life_table, number_of_years_to_shift, file_name_to_print)

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
survivorship_curve_vector = calc_survivorship_curve_vector(count_of_age_cohorts, survival_vector);
raw_growth_rate = sum(original_fecundity_vector .* survivorship_curve_vector);

fecundity_func = fit(age_vector,original_fecundity_vector,'linearinterp') ;
startValue = 1; % age 1
endValue = count_of_age_cohorts;
new_length_of_adulthood = count_of_age_cohorts - number_of_years_to_shift;
% new sampling points for interpolating new sample fecundties
new_adult_fecundity_vector_sample_points = linspace(startValue, endValue, new_length_of_adulthood);
% interpolated values for fecundities of the new adult age classes
new_adult_fecundity_vector_values = fecundity_func(new_adult_fecundity_vector_sample_points);

new_fecundity_vector = [zeros(number_of_years_to_shift,1); new_adult_fecundity_vector_values];

recalculated_growth_rate = sum(new_fecundity_vector .* survivorship_curve_vector);

% scale up new fecundity vector so new growth rate is same as the old growth rate
new_fecundity_vector = new_fecundity_vector .* (raw_growth_rate/recalculated_growth_rate);


life_table(:,3) = new_fecundity_vector;

col_header={'age','s','b'}; 

xlswrite(file_name_to_print, col_header, 'Sheet1', 'A1');
xlswrite(file_name_to_print, life_table, 'Sheet1', 'A2');

end

% This formula produces the l-sub-x vector descrbied in Waples 2011)
function [survivorship_curve_vector] = calc_survivorship_curve_vector(count_of_age_cohorts, survival_vector);

survivorship_curve_vector = ones(count_of_age_cohorts, 1);

for i = 2:count_of_age_cohorts
    survivorship_curve_vector(i) = survivorship_curve_vector(i - 1) * survival_vector(i - 1);
end
end