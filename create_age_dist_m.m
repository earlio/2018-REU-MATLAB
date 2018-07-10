function [age_dist_m] = create_age_dist_m(number_generations, population_0, leslie_matrix)

time = 1:number_generations; %creates a time vector using the set number of generations

age_dist_m = zeros(length(population_0), length(time)); %creates a matrix with age distributions at each time step
age_dist_m(:,1) = population_0; %sets the initial age distributions to the initial population values in population_0

total_population_0 = sum(population_0); %calculates the total population at time = 0 by adding the values of population_0
total_population_v = zeros(1, length(time)); %initializes a vector of total population values
total_population_v(1) = total_population_0; %sets the first entry of the total population vector equal to the initial total population

for i = 2:length(time)
    age_dist_m(:,i) = round(leslie_matrix*age_dist_m(:,i-1)); %applies the leslie matrix to the previous age distribution for each time step
    
    total_population_v(i) = sum(age_dist_m(:,i)); %adds the total population at time step i to the total population vector
    if isequal(
end