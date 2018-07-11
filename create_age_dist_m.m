function [age_dist_m] = create_age_dist_m(number_generations, population_0, leslie_matrix)
time = 1:number_generations; %creates a time vector using the set number of generations

total_population_0 = sum(population_0); %calculates the total population at time = 0 by adding the values of population_0
total_population_v = zeros(1, length(time)); %initializes a vector of total population values
total_population_v(1) = total_population_0; %sets the first entry of the total population vector equal to the initial total population

%% Burn-In %%

burn_in_m = zeros(length(population_0),50);
for i = 2:50
    burn_in_m(:,i) = round(leslie_matrix*burn_in_m(:,i-1));
    total_population_v(i) = sum(burn_in_m(:,i));
    
    if isequal(burn_in_m(:,i),burn_in_m(:,i-1))
        age_dist_m = burn_in_m;
        burn_in = true;
        break
    else
        burn_in = false;
    end
    if total_population_v(i) >= 10^10
        age_dist_m = burn_in_m(:,1:i);
        disp("Warning: Population size is too high") %ends the creation of the matrix
        burn_in = true;
        break
    end
end



%% Matrix Creation %%
if burn_in == false

age_dist_m = zeros(length(population_0), length(time)); %creates a matrix with age distributions at each time step
age_dist_m(:,1) = burn_in_m(:,end); %sets the initial age distributions to the initial population values in population_0

for i = 2:length(time)
    age_dist_m(:,i) = round(leslie_matrix*age_dist_m(:,i-1)); %applies the leslie matrix to the previous age distribution for each time step
    
    total_population_v(i) = sum(age_dist_m(:,i)); %adds the total population at time step i to the total population vector
    
    if total_population_v(i) >= 10^10
        age_dist_m = age_dist_m(:,1:i);
        disp("Warning: Population size is too high") %ends the creation of the matrix
        break
    end
    if isequal(age_dist_m(:,i),age_dist_m(:,i-1))
        if isequal(age_dist_m(:,i),age_dist_m(:,i-2)) && isequal(age_dist_m(:,i),age_dist_m(:,i-3))
        age_dist_m = age_dist_m(:,1:i);
        break %breaks the loop at a steady state
        end
    end
end
end