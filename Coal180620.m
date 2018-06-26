%% Short Hand %%

% "_0" label is time 0
% "_m" label indicates a matrix
% "_v" label indicates a vector
% "dist" = distribution
%% User Inputs %%

% population_0 = input('Input an initial population vector and press enter: ');
% number_generations = input('Input a number of generations and press enter: ');
% leslie_matrix = input('Input a Leslie Matrix and press enter: ');
% lineage_count = input('Input a number of lineages to track and press enter: ');

%% Example Inputs %% 

%EXAMPLE 1
% population_0 = [2, 2, 2]; %example initial population vector
% number_generations = 10; %example number of generations
% leslie_matrix = [0 1 1.1; 0.6 0 0; 0 0.5 0]; %example leslie matrix
% lineage_count = 2; %example number of lineages to track

%EXAMPLE 2
population_0 = [2, 2, 2, 5]; %example initial population vector
number_generations = 20; %example number of generations
leslie_matrix = [0 1 1.1 1.2; 0.6 0 0 0; 0 0.5 0 0; 0 0 0.25 0]; %sample leslie matrix
lineage_count = 2; %example number of lineages to track 

%EXAMPLE 3
% population_0 = [56 72 16 84 23];
% number_generations = 100; 
% leslie_matrix = [0 0.5 0.75 1 1.25; 0.9 0 0 0 0; 0 0.75 0 0 0; 0 0 0.6 0 0; 0 0 0 0.4 0];
% lineage_count = 2;

%EXAMPLE 4
% population_0 = [506 823 234 348 294 1844]
% number_generations = 50;
% leslie_matrix = [0 0.1 0.1 0.1 0.1 0; 0.95 0 0 0 0 0; 0 0.92 0 0 0 0; 0 0 0.87 0 0 0; 0 0 0 0.56 0 0; 0 0 0 0 0.3 0];
% lineage_count = 2; 


%% Establish Age Distribution Matrix (age_dist_m) %%

time = 1:number_generations; %creates a time vector using the set number of generations

age_dist_m = zeros(length(population_0), length(time)); %creates a matrix with age distributions at each time step
age_dist_m(:,1) = population_0; %sets the initial age distributions to the initial population values in population_0

total_population_0 = sum(population_0); %calculates the total population at time = 0 by adding the values of population_0
total_population_v = zeros(1, length(t)); %initializes a vector of total population values
total_population_v(1) = total_population_0; %sets the first entry of the total population vector equal to the initial total population

for i = 2:length(t)
    age_dist_m(:,i) = round(leslie_matrix*age_dist_m(:,i-1)); %applies the leslie matrix to the previous age distribution for each time step
    total_population_v(i) = sum(age_dist_m(:,i)); %adds the total population at time step i to the total population vector
end

%% Establish a Matrix of Individuals (individuals_m) %%

max_population = max(total_population_v); %sets the column dimension of the individuals matrix by determining the maximum population value from the total population vector

individuals_m = -1*ones(number_generations,max_population); %creates a matrix with number_generations rows and max_population columns where every entry is -1

for i = 1:length(t)
    individuals_m(i,1:age_dist_m(1,i)) = 0; %boundary case, set the first n individuals to be age zero (using the first index of the age_dist_m matrix)
    for j = 2:length(population_0)
        index = 1+sum(age_dist_m(1:j-1,i)); %determine the number of individuals already assigned to an age, "+1" (MATLAB indices are inclusive)
        individuals_m(i,index:(index-1+age_dist_m(j,i))) = j-1; %assign individuals of the next age group to the newest open spaces in the AgeM matrix
    end
end


%% Extract the Terminal Population and Choose Number of Individuals to Track  %%

terminal_population = extract_terminal_population(age_dist_m); %returns the portion of the last row which contains individuals
lineage_count = 2; %number of lineages user wants to track, will be updated with user input later. 
indices = terminal_indices(terminal_population, lineage_count); %function which returns two random indices from the final row of individuals. 

%% Track the lineages to an MRCA %%


    
    
    


%  
% disp(terminal_population);
% % ez we'll start with just two hard-coded pairs to find the mrca of.
% 
% lineage_a_current_age = 0;
% lineage_b_current_age = 0;
% 
% 
% mrca = calculate_mrca(lineage_a_current_age, lineage_b_current_age, ...
%     terminal_population, life_table, numgen);
% 
% disp(mrca)
% 
%     
%     
% 
