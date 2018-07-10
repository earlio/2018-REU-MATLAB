%% top level %%
%this file organizes a test of the generation-by-generation coalescent.
%Input a test life table excel file, an initial population vector, a number
%of generations, a lineage count, and number of
%iterations to generate a histogram of results for t-mrca

%% life_to_leslie calls %%

%sample_choice = input('Input a number to choose a sample population. 1 = Sheep, 2 = Mice, 3 = Sardines, 4 = Primroses');
sample_choice = 4; %chooses which sample to test, change this value to change the sample, 1 = Sheep, 2 = Mice, 3 = Sardines, 4 = Primroses

if isequal(sample_choice,1)
[leslie_matrix, ages, alpha1_AL_lambda]= life_to_leslie('Sample_LT1.xlsx','A3:C15');

elseif isequal(sample_choice,2)
[leslie_matrix, ages, alpha1_AL_lambda] = life_to_leslie('Sample_LT2.xlsx','A3:C8');

elseif isequal(sample_choice,3)
[leslie_matrix, ages, alpha1_AL_lambda] = life_to_leslie('Sample_LT3.xlsx','A3:C15');

elseif isequal(sample_choice,4)
[leslie_matrix, ages, alpha1_AL_lambda] = life_to_leslie('Sample_LT4.xlsx','A3:C24');

end

%% Normalize the Leslie Matrix %%

% lambda = eig(leslie_matrix);
% 
% if lambda(1)> 0
%     scaling = lambda(1);
% else
%     scaling = lambda(2);
% end
% 
% leslie_matrix = leslie_matrix./scaling;



%% Set Initial Parameters %%

% population_0 = input('Input an initial population vector and press enter: ');
% number_generations = input('Input a number of generations and press enter: ');
% lineage_count = input('Input a number of lineages to track and press enter: ');

% number_iterations = input('Input a number of iterations and press enter: ');


population_0 = 100*ones(1,ages); %initial population vector, sets each initial population to 10 for now.
number_generations = 200; %number of generations,
lineage_count = 2; %number of lineages to track, k=2 for now
number_iterations = 1; %number of times the functions are run

number_iterations_v = zeros(2,number_iterations); %creates a vector to keep track of data for each iteration

%% Create the demographic matrix %%

age_dist_m = create_age_dist_m(number_generations, population_0, leslie_matrix); %%calls the function to create the demographic matrix

%% Choose starting indices and establish the genealogy matrix %%



for j = 1:2
if isequal(j,1)
    age_i = -1; %testing for Ne
elseif isequal(j,2)
    age_i = 0; %testing for Nb
end
for n = 1:number_iterations
initial_values = terminal_indices(lineage_count, age_dist_m, age_i); %function which returns two random indices from the final row of individuals. 

genealogy_m = -1*ones(number_generations, lineage_count, 2); %initialize the 3-D genealogy matrix
genealogy_m(end,:,1) = initial_values(1,:); genealogy_m(end,:,2) = initial_values(2,:); %set the front row to the indices and the back row to the ages specified in the initial_values matrix

[mrca, complete_genealogy,coal_events] = calc_mrca_b(genealogy_m, leslie_matrix, age_dist_m)

if isequal(age_i,-1)
    number_iterations_v(1,n) = mrca; %sets the first row of number iterations to Ne results
elseif isequal(age_i,0)
    number_iterations_v(2,n) = mrca; %sets second row of number iterations to Nb results
end
end
end

number_iterations_v