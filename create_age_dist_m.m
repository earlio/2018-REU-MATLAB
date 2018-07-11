function [age_dist_m] = create_age_dist_m(number_generations, population_0, leslie_matrix,burn_in_gens)
time = 1:number_generations; %creates a time vector using the set number of generations

total_population_0 = sum(population_0); %calculates the total population at time = 0 by adding the values of population_0
total_population_v = zeros(1, length(time)); %initializes a vector of total population values
total_population_v(1) = total_population_0; %sets the first entry of the total population vector equal to the initial total population

%% Burn-In %%
warning = 0;
burn_in_m = zeros(length(population_0),50);
burn_in_m(:,1) = population_0;
for i = 2:burn_in_gens
    burn_in_m(:,i) = round(leslie_matrix*burn_in_m(:,i-1));
    total_population_v(i) = sum(burn_in_m(:,i));
    
    if isequal(burn_in_m(:,i),burn_in_m(:,i-1)) && isequal(warning,0)
        fprintf("The Population reaches a steady state at burn-in t = %d \n", i);
        burn_in = false;
        warning = 1;
    else
        burn_in = false;
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

for i = 2:length(time)
    age_dist_m(:,i) = round(leslie_matrix*age_dist_m(:,i-1)); %applies the leslie matrix to the previous age distribution for each time step
    
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