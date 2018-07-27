function [age_dist_m, options_m] = create_age_dist_m(number_generations, population_0, leslie_matrix,burn_in_gens)
%Inputs:
%1. number_generations - The preset number of generations that we are testing the data on
%2. population_0 - a vector with the initial population of individuals for
%each age class
%3. leslie_matrix - a leslie matrix corresponding to the input life table,
%generated using the life_to_leslie function
%4. burn-in-gens - the number of generations required for the population to
%stabilize

% This function generates a burn-in matrix to allow the population to
% stabilize, and an age distribution matrix which contains the distribution
% of individuals among the age classes over the preset number of
% generations. The age distribution matrix begins with the final column of
% the burn-in matrix. This function also displays warnings if the
% population size of an age class grows too large, and a warning of when
% the population stabilizes to ensure the number of burn-in generations is
% correctly set.


%Outputs:
%1. age_dist_m - The only output is the demographic matrix/age distribution matrix



time = 1:number_generations; %creates a time vector using the set number of generations

total_population_0 = sum(population_0); %calculates the total population at time = 0 by adding the values of population_0
total_population_v = zeros(1, length(time)); %initializes a vector of total population values
total_population_v(1) = total_population_0; %sets the first entry of the total population vector equal to the initial total population

%% Burn-In %%
warning = 0;
burn_in_m = zeros(length(population_0),burn_in_gens);
burn_in_m(:,1) = population_0;
size(burn_in_m)
size(leslie_matrix)
for i = 2:burn_in_gens
    burn_in_m(:,i) = round(leslie_matrix*burn_in_m(:,i-1));
    total_population_v(i) = sum(burn_in_m(:,i));
    
    if isequal(burn_in_m(:,i),burn_in_m(:,i-1)) && isequal(warning,0)
        fprintf("The Population reaches a steady state at burn-in t = %d \n", i);
        burn_in = false;
        warning = 1;
    else
        burn_in = false;
        fprintf("The Population has not reached steady state after %d burn-in generations \n", burn_in_gens);
    end
    size_check_v = burn_in_m(:,i) > (2^50)*ones(size(burn_in_m,1),1);
    if isequal(size_check_v,zeros(size(burn_in_m,1),1) == 0) && (isequal(warning,0))
        fprintf("The Population is too large at burn-in t = %d \n",i);
        burn_in = false;
        warning = 1;
    end
end

%% Matrix Creation %%
if burn_in == false
    
    age_dist_m = zeros(length(population_0), length(time)); %creates a matrix with age distributions at each time step
    age_dist_m(:,1) = burn_in_m(:,end); %sets the initial age distributions to the initial population values in population_0
    
    % ez options_m is a lookup table for the cumululative age ranges
    % used to assign parents in calc_mrca_b and should speed execution by a
    % lot. options_m(:,;,1) is "lower" and options_m(:,:,2) is "upper".

    age_cohort_count = length(population_0);
    options_m = zeros(age_cohort_count, length(time), 2);
    options_m(:,1,:) =  build_options_vectors(age_dist_m(:,1), age_cohort_count);
    
    for i = 2:length(time)
    
%       age_dist_m(:,i) = round(leslie_matrix*age_dist_m(:,i-1)); %applies the leslie matrix to the previous age distribution for each time step

        rounded_demo_vector = round(leslie_matrix*age_dist_m(:,i-1)); %applies the leslie matrix to the previous age distribution for each time step
        age_dist_m(:,i) = rounded_demo_vector;
        options_m(:,i,:) =  build_options_vectors(rounded_demo_vector, age_cohort_count);
        
        total_population_v(i) = sum(age_dist_m(:,i)); %adds the total population at time step i to the total population vector
        size_check_v = age_dist_m(:,i) > (2^50)*ones(size(age_dist_m,1),1);
        if isequal(size_check_v,zeros(size(age_dist_m,1),1) == 0) && (isequal(warning,0))
            fprintf("The Population is too large at t = %d \n",i);
            warning = 1;
        end
        if isequal(age_dist_m(:,i),age_dist_m(:,i-1)) && isequal(warning,0)
            age_dist_m = age_dist_m(:,1:i);
            fprintf("The Population reaches a steady state at t = %d \n", i);
            warning = 1;
        end
    end
end
end

function annual_options_vector = build_options_vectors(age_dist_slice, age_cohort_count)

annual_options_vector = zeros(age_cohort_count,1,2);
lower_subtotal = 1;
upper_subtotal = age_dist_slice(1,1);
annual_options_vector(1, 1, 1) = lower_subtotal;
annual_options_vector(1, 1, 2) = upper_subtotal;

for i = 2:age_cohort_count
    %lower
    lower_subtotal = lower_subtotal + age_dist_slice(i - 1);
    annual_options_vector(i, 1, 1) = lower_subtotal;
    %upper
    upper_subtotal = upper_subtotal + age_dist_slice(i);
    annual_options_vector(i, 1, 2) = upper_subtotal;
    
end

end