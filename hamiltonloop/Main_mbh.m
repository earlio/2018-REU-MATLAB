

rng('shuffle'); % random seed for random number generator

% declare main glaobal variables for program
total_pop_N = 500; % size of population for all age classes

number_generations = 5000; % number of generations

burn_in_gens = 50; % number of generations for burn in of population growth

lineage_count = 2; % number of lineages to sample to determine time to MRCA 

iterations = 500; % number of iterations of sampling from population


fprintf('----------------------------------------------------\n');
fprintf('Simulation of time to MRCA in an age-structured coalescent\n\n');

t = datetime('now','TimeZone','local','Format','d-MMM-y HH:mm:ss Z'); % get date and time
date_string = datestr(t); % convert date and time to string
fprintf('%s\n\n',date_string);


%% open file with life table, get Leslie matrix for population growth

% provide path name to life table file
file_path_name = 'Sample_LT1.xlsx';
fprintf('Life table file: %s\n\n',file_path_name);

% rescale life table so that population size is constant through time
rescale = true;

[leslie_matrix,age_classes,scaling, life_table_m] = life_to_leslie(file_path_name,rescale);

if rescale == true
    fprintf('original life table had growth rate of %f\n', scaling);
    fprintf('life table has been rescaled to constant population size over time\n\n');
else
    fprintf('life table has a growth rate of %f\n\n', scaling);
end

%% Establish Age Distribution Matrix (age_dist_m)

time = 1:number_generations; % create a time vector using the number of generations

uniform_age_cohort_N = round(total_pop_N/(age_classes - 1)); % initial uniform size for all age cohorts

age_dist_m = zeros(age_classes, number_generations); % initialize a matrix with age-specific N at each time step
temp_age_dist_m = zeros(age_classes, burn_in_gens); % initialize a matrix with age-specific N at each time step

temp_age_dist_m(1:(age_classes-1),1) = uniform_age_cohort_N; % sets the initial age cohort sizes to be uniform

total_population_0 = sum(temp_age_dist_m(:,1)); % compute the total population at time = 0 by summing all elements in the first column of age_dist_m

% check to see if all age cohorts sum to total population size,
% if not make adjustment
if total_population_0 - total_pop_N > 0
    
    temp_age_dist_m(1,1) = temp_age_dist_m(1,1) + (total_population_0 - total_pop_N); 
    
elseif total_population_0 - total_pop_N < 0
    
     temp_age_dist_m(1,1) =  temp_age_dist_m(1,1) + (total_pop_N - total_population_0);
     
end

% burn in by carrying out population growth for a few generations to remove initilization effects
for i = 2:burn_in_gens
    temp_age_dist_m(:,i) = round(leslie_matrix*temp_age_dist_m(:,i-1)); % applies the leslie matrix to the previous age distribution for each time step
end

% start age distribution matrix with final generation of burn in temporary
% population
age_dist_m(:,1) = temp_age_dist_m(:,burn_in_gens); 


total_population_N = zeros(1, length(time)); %initialize a vector of total population size at each time

total_population_N(1,1) = total_pop_N; % set the first entry of the total population size vector equal to the initial total population size

for i = 2:number_generations
    
    age_dist_m(:,i) = round(leslie_matrix*age_dist_m(:,i-1)); % applies the leslie matrix to the previous age distribution for each time step
    
    total_population_N(i) = sum(age_dist_m(:,i)); % adds the total population at time step i to the total population vector
end


% %% Extract the Terminal Population and Choose Number of Lineages to Track
% 
% %terminal_population = extract_terminal_population2(age_dist_m); %returns the portion of the last row which contains individuals
% 
% % an integer within the range of age classes to sample individuals in that age class. 
% % Enter -1 to sample individuals of randon ages. 
% age_i = -1; 
% 
% initial_values = terminal_indices(lineage_count,age_dist_m,age_i); % function which samples lineages from all lineages in the present. 
% 
% genealogy_m = -1*ones(number_generations, lineage_count, 2); % initialize the 3-D genealogy matrix *** need to describe the rows, cols and pages!
% 
% genealogy_m(end,:,1) = initial_values(1,:); genealogy_m(end,:,2) = initial_values(2,:); % set the front row to the indices and the back row to the ages specified in the initial_values matrix
% 
% disp(genealogy_m);
% %% Track the lineages to an MRCA
% 
% [mrca,complete_genealogy,coal_events] = calc_mrca_b(genealogy_m, leslie_matrix, age_dist_m);
% 
% if isequal(mrca, number_generations)
%     disp('no mrca')
% else
%     mrca
% end
% %complete_genealogy
% %coal_events


% sample pairs of random age lineages
age_i = -1; % set lineages sampled in present to random ages

no_mrca_random = 0; % counter for number of times no MRCA was found

mrca_random = zeros(1,iterations); % allocate space for results

for iter=1:iterations

    initial_values = terminal_indices(lineage_count,age_dist_m,age_i); % samples lineage from all lineages in the present. 

    genealogy_m = -1*ones(number_generations, lineage_count, 2); % initialize the 3-D genealogy matrix *** need to describe the rows, cols and pages!

    genealogy_m(end,:,1) = initial_values(1,:); genealogy_m(end,:,2) = initial_values(2,:); % set the front row to the indices and the back row to the ages specified in the initial_values matrix

    % Track the lineages to an MRCA

    [mrca,complete_genealogy,coal_events] = calc_mrca_ez(genealogy_m, life_table_m, age_dist_m);
 %   [mrca,complete_genealogy,coal_events] = calc_mrca_b(genealogy_m, leslie_matrix, age_dist_m);

    if mrca == number_generations
        no_mrca_random = no_mrca_random + 1; % increment counter
    else
        mrca_random(1,iter) = mrca;
        
    end
end % for iter

    % compute mean time to MRCA without iterations that did not experience
    % coalescence
    num_non_zero_elements = iterations - no_mrca_random;
    sum_mrca_random = sum(mrca_random);
    mean_random = sum_mrca_random/num_non_zero_elements;
    median_random = median(mrca_random(1,1:num_non_zero_elements));


% sample pairs of age zero lineages
age_i = 0; % set lineages sampled in present to age zero

no_mrca_zero = 0; % counter for number of times no MRCA was found

mrca_zero = zeros(1,iterations); % allocate space for results

for iter=1:iterations

    iter
    
    initial_values = terminal_indices(lineage_count,age_dist_m,age_i); % samples lineage from all lineages in the present. 

    genealogy_m = -1*ones(number_generations, lineage_count, 2); % initialize the 3-D genealogy matrix *** need to describe the rows, cols and pages!

    genealogy_m(end,:,1) = initial_values(1,:); genealogy_m(end,:,2) = initial_values(2,:); % set the front row to the indices and the back row to the ages specified in the initial_values matrix

    % Track the lineages to an MRCA

    [mrca,complete_genealogy,coal_events] = calc_mrca_ez(genealogy_m, life_table_m, age_dist_m);
 %   [mrca,complete_genealogy,coal_events] = calc_mrca_b(genealogy_m, leslie_matrix, age_dist_m);

    if mrca == number_generations
        no_mrca_zero = no_mrca_zero + 1; % increment counter
    else
        mrca_zero(1,iter) = mrca;
    end
end % for iter

    % compute mean time to MRCA without iterations that did not experience
    % coalescence
    num_non_zero_elements = iterations - no_mrca_zero;
    sum_mrca_zero = sum(mrca_zero);
    mean_zero = sum_mrca_zero/num_non_zero_elements;
    median_zero = median(mrca_zero(1,1:num_non_zero_elements));


    figure;

    hold on;
    subplot(2,1,1);
    hist(mrca_random); 
    xlabel('time to MRCA - random lineage pairs')
    ylabel('Count')

    subplot(2,1,2);
    hist(mrca_zero); 
    xlabel('time to MRCA - age zero lineage pairs')
    ylabel('Count')

    suptitle('Distributions of coalescence times');
    hold off;

    
    
    figure;
    hold on;
    subplot(1,2,1);
    boxplot(mrca_zero, 'Labels',{'age zero lineage pairs'});
    
    subplot(1,2,2);
    boxplot(mrca_random, 'Labels',{'random age lineage pairs'});
    hold off;

    
    fprintf('Summary of simulation:\n\n');
        
    fprintf('Number of iterations was %i\n', iterations);
    fprintf('Maximum number of generations was %i\n', number_generations);
    fprintf('Total population size was %i\n', total_pop_N);
    
    fprintf('\n');
    
    fprintf('zero age lineages coalescence times:\n');
    fprintf('average: %f \n', mean_zero);
    fprintf('median: %f \n', median_zero);
    
    fprintf('\n');
   
    fprintf('random age lineages coalescence times:\n');
    fprintf('average: %f \n', mean_random);
    fprintf('median: %f \n', median_random);    
    




fprintf('----------------------------------------------------\n');



function [parent] = sample_lineage(age)


    % sample one lineage given an age class and a generation
    age = parent_age_class; % age of lineage to find
    
    generation = numgen - 1; % generation of lineage to find
    
    found_one = 0;
    while found_one == 0

        index = randi(total_pop_sizes(generation,1)); % get random integer 1:total population size at time numgen
        
        if pop_in_present(1,index) == age % lineage must be given age 
            parent = index;
            found_one = 1;
        else
            found_one = 0; % keep looking for age zero individual
        end
    end % while

    disp(parent);
    
end

