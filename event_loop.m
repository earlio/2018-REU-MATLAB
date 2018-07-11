%% top level %%
%this file organizes a test of the generation-by-generation coalescent.
%Input a test life table excel file, an initial population vector, a number
%of generations, a lineage count, and number of
%iterations to generate a histogram of results for t-mrca




%% Set Up Output %%


fprintf('----------------------------------------------------\n');
fprintf('Simulation of time to MRCA in an age-structured coalescent\n\n');

t = datetime('now','TimeZone','local','Format','d-MMM-y HH:mm:ss Z'); % get date and time
date_string = datestr(t); % convert date and time to string
fprintf('%s\n\n',date_string);

%% Set Initial Parameters %%

% population_0 = input('Input an initial population vector and press enter: ');
% number_generations = input('Input a number of generations and press enter: ');
% lineage_count = input('Input a number of lineages to track and press enter: ');

% number_iterations = input('Input a number of iterations and press enter: ');


population_0 = 100*ones(1,ages); %initial population vector, sets each initial population to 10 for now.
number_generations = 300; %number of generations,
lineage_count = 2; %number of lineages to track, k=2 for now
number_iterations = 1; %number of times the functions are run
burn_in_gens  = 50; 


number_iterations_v = zeros(2,number_iterations); %creates a vector to keep track of data for each iteration


%% Leslie Matrix Creation and Normalization %%

% provide path name to life table file
file_path_name = 'WF_life_table.xlsx';
fprintf('Life table file: %s\n\n',file_path_name);

% rescale life table so that population size is constant through time
rescale = false;

[leslie_matrix,age_classes,scaling] = life_to_leslie_2(file_path_name,rescale);

if rescale == true
    fprintf('original life table had growth rate of %f\n', scaling);
    fprintf('life table has been rescaled to constant population size over time\n\n');
else
    fprintf('life table has a growth rate of %f\n\n', scaling);
end


%% Create the demographic matrix %%

age_dist_m = create_age_dist_m(number_generations, population_0, leslie_matrix, burn_in_gens); %%calls the function to create the demographic matrix

%% Choose starting indices and establish the genealogy matrix %%


for j = 1:2
if isequal(j,1)
    age_i = -1; %testing for Ne
elseif isequal(j,2)
    age_i = 0; %testing for Nb
end
for n = 1:number_iterations
initial_values = terminal_indices(lineage_count, age_dist_m, age_i) %function which returns two random indices from the final row of individuals. 

genealogy_m = -1*ones(number_generations, lineage_count, 2); %initialize the 3-D genealogy matrix
genealogy_m(end,:,1) = initial_values(1,:); genealogy_m(end,:,2) = initial_values(2,:); %set the front row to the indices and the back row to the ages specified in the initial_values matrix

[mrca, complete_genealogy,coal_events] = calc_mrca_b(genealogy_m, leslie_matrix, age_dist_m);

if isequal(age_i,-1)
    number_iterations_v(1,n) = mrca; %sets the first row of number iterations to Ne results
elseif isequal(age_i,0)
    number_iterations_v(2,n) = mrca; %sets second row of number iterations to Nb results
end
end
end

number_iterations_v