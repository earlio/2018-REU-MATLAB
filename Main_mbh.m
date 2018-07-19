

rng('shuffle'); % random seed for random number generator

% declare main glaobal variables for program
total_pop_N = 500; % size of population for all age classes

<<<<<<< HEAD
number_generations = 5000; % number of generations

burn_in_gens = 50; % number of generations for burn in of population growth
=======
number_generations = 3000; % number of generations

burn_in_gens = 102; % number of generations for burn in of population growth
>>>>>>> L2Lv5

lineage_count = 2; % number of lineages to sample to determine time to MRCA 

iterations = 500; % number of iterations of sampling from population

<<<<<<< HEAD

=======
>>>>>>> L2Lv5
fprintf('----------------------------------------------------\n');
fprintf('Simulation of time to MRCA in an age-structured coalescent\n\n');

t = datetime('now','TimeZone','local','Format','d-MMM-y HH:mm:ss Z'); % get date and time
date_string = datestr(t); % convert date and time to string
fprintf('%s\n\n',date_string);

<<<<<<< HEAD
=======
fprintf('Initial total population size: %g\n\n', total_pop_N);
>>>>>>> L2Lv5

%% open file with life table, get Leslie matrix for population growth

% provide path name to life table file
<<<<<<< HEAD
file_path_name = 'Sample_LT4.xlsx';
fprintf('Life table file: %s\n\n',file_path_name);

% rescale life table so that population size is constant through time
rescale = false;

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
=======
file_path_name = 'Sample_sage_grouse_life_table.xlsx';
fprintf('Life table file: %s\n\n',file_path_name);

% scaling factor to adjust population growth rate; 
% set = 1/lambda for same lambda as in the life table
% when scale = 1.0 there is no population growth (lambda = 1.0)
% when scale > 1.0 there is negative population growth rate (lambda < 1.0)
% when scale < 1.0 there is positive population growth rate (lambda > 1.0)

%scale = 0.98;
scale = 1.0;

[life_table_m,leslie_matrix,age_classes,orig_lambda,mod_lambda,CV_fecundity,G,alpha,AL] = life_to_leslie(file_path_name,scale);

fprintf('life table has a population growth rate of %f\n', orig_lambda);
fprintf('(set scale to 1/lambda or = %f for original rate)\n\n', 1/orig_lambda);

fprintf('scaling factor for population growth rate was %f\n', scale);
fprintf('after rescaling, the population growth rate is %f\n\n', mod_lambda);

fprintf('number of age classes (maximum life span, omega): %g\n', age_classes);
fprintf('G, average age of a parent: %g\n', G);
fprintf('alpha, age of first reproductuve maturity: %g\n\n', alpha);

fprintf('coefficient of variation in age-specific fecundity: %f\n\n', CV_fecundity);
%% Establish Age Distribution Matrix (age_dist_m)

time = 1:number_generations; % create a time vector using the number of generations
 
uniform_age_cohort_N = round(total_pop_N/(age_classes)); % initial uniform size for all age cohorts
 
population_0 = zeros(1,age_classes);
 
population_0(1:(age_classes)) = uniform_age_cohort_N;
 
total_population_0 = sum(population_0(:,1));
 
 
if total_population_0 - total_pop_N > 0
    
    population_0(1) = population_0(1) + (total_population_0 - total_pop_N); 
    
elseif total_population_0 - total_pop_N < 0
    
     population_0(1) =  population_0(1) + (total_pop_N - total_population_0);
     
end

% create the demographic matrix of population size for each age class over time
age_dist_m = create_age_dist_m(number_generations, population_0, leslie_matrix, burn_in_gens); 
>>>>>>> L2Lv5


% sample pairs of random age lineages
age_i = -1; % set lineages sampled in present to random ages

no_mrca_random = 0; % counter for number of times no MRCA was found

mrca_random = zeros(1,iterations); % allocate space for results

for iter=1:iterations

    initial_values = terminal_indices(lineage_count,age_dist_m,age_i); % samples lineage from all lineages in the present. 

    genealogy_m = -1*ones(number_generations, lineage_count, 2); % initialize the 3-D genealogy matrix *** need to describe the rows, cols and pages!

    genealogy_m(end,:,1) = initial_values(1,:); genealogy_m(end,:,2) = initial_values(2,:); % set the front row to the indices and the back row to the ages specified in the initial_values matrix

    % Track the lineages to an MRCA

<<<<<<< HEAD
    [mrca,complete_genealogy,coal_events] = calc_mrca_ez(genealogy_m, life_table_m, age_dist_m);
 %   [mrca,complete_genealogy,coal_events] = calc_mrca_b(genealogy_m, leslie_matrix, age_dist_m);
=======
    [mrca,complete_genealogy,coal_events] = calc_mrca_b(genealogy_m, leslie_matrix, age_dist_m);
>>>>>>> L2Lv5

    if mrca == number_generations
        no_mrca_random = no_mrca_random + 1; % increment counter
    else
        mrca_random(1,iter) = mrca;
        
    end
end % for iter

    % compute mean time to MRCA without iterations that did not experience
    % coalescence
<<<<<<< HEAD
    num_non_zero_elements = iterations - no_mrca_random;
    sum_mrca_random = sum(mrca_random);
    mean_random = sum_mrca_random/num_non_zero_elements;
    median_random = median(mrca_random(1,1:num_non_zero_elements));
=======
    num_non_zero_sims_random = iterations - no_mrca_random;
    sum_mrca_random = sum(mrca_random);
    mean_random = sum_mrca_random/num_non_zero_sims_random;
    median_random = median(mrca_random(1,1:num_non_zero_sims_random));
>>>>>>> L2Lv5


% sample pairs of age zero lineages
age_i = 0; % set lineages sampled in present to age zero

no_mrca_zero = 0; % counter for number of times no MRCA was found

mrca_zero = zeros(1,iterations); % allocate space for results

for iter=1:iterations

<<<<<<< HEAD
    iter
    
=======
>>>>>>> L2Lv5
    initial_values = terminal_indices(lineage_count,age_dist_m,age_i); % samples lineage from all lineages in the present. 

    genealogy_m = -1*ones(number_generations, lineage_count, 2); % initialize the 3-D genealogy matrix *** need to describe the rows, cols and pages!

<<<<<<< HEAD
    genealogy_m(end,:,1) = initial_values(1,:); genealogy_m(end,:,2) = initial_values(2,:); % set the front row to the indices and the back row to the ages specified in the initial_values matrix

    % Track the lineages to an MRCA

    [mrca,complete_genealogy,coal_events] = calc_mrca_ez(genealogy_m, life_table_m, age_dist_m);
 %   [mrca,complete_genealogy,coal_events] = calc_mrca_b(genealogy_m, leslie_matrix, age_dist_m);
=======
    genealogy_m(end,:,1) = initial_values(1,:);
    genealogy_m(end,:,2) = initial_values(2,:); % set the front row to the indices and the back row to the ages specified in the initial_values matrix

    % Track the lineages to an MRCA

    [mrca,complete_genealogy,coal_events] = calc_mrca_b(genealogy_m, leslie_matrix, age_dist_m);
>>>>>>> L2Lv5

    if mrca == number_generations
        no_mrca_zero = no_mrca_zero + 1; % increment counter
    else
        mrca_zero(1,iter) = mrca;
    end
end % for iter

    % compute mean time to MRCA without iterations that did not experience
    % coalescence
<<<<<<< HEAD
    num_non_zero_elements = iterations - no_mrca_zero;
    sum_mrca_zero = sum(mrca_zero);
    mean_zero = sum_mrca_zero/num_non_zero_elements;
    median_zero = median(mrca_zero(1,1:num_non_zero_elements));
=======
    num_non_zero_sims_zero = iterations - no_mrca_zero;
    sum_mrca_zero = sum(mrca_zero);
    mean_zero = sum_mrca_zero/num_non_zero_sims_zero;
    median_zero = median(mrca_zero(1,1:num_non_zero_sims_zero));
>>>>>>> L2Lv5


    figure;

    hold on;
    subplot(2,1,1);
<<<<<<< HEAD
    hist(mrca_random); 
=======
    hist(mrca_random(1,1:num_non_zero_sims_random)); 
>>>>>>> L2Lv5
    xlabel('time to MRCA - random lineage pairs')
    ylabel('Count')

    subplot(2,1,2);
<<<<<<< HEAD
    hist(mrca_zero); 
=======
    hist(mrca_zero(1,1:num_non_zero_sims_zero)); 
>>>>>>> L2Lv5
    xlabel('time to MRCA - age zero lineage pairs')
    ylabel('Count')

    suptitle('Distributions of coalescence times');
    hold off;

    
    
    figure;
    hold on;
    subplot(1,2,1);
<<<<<<< HEAD
    boxplot(mrca_zero, 'Labels',{'age zero lineage pairs'});
    
    subplot(1,2,2);
    boxplot(mrca_random, 'Labels',{'random age lineage pairs'});
=======
    boxplot(mrca_zero(1,1:num_non_zero_sims_zero), 'Labels',{'age zero lineage pairs'});
    
    subplot(1,2,2);
    boxplot(mrca_random(1,1:num_non_zero_sims_random), 'Labels',{'random age lineage pairs'});
>>>>>>> L2Lv5
    hold off;

    
    fprintf('Summary of simulation:\n\n');
        
    fprintf('Number of iterations was %i\n', iterations);
    fprintf('Maximum number of generations was %i\n', number_generations);
    fprintf('Total population size was %i\n', total_pop_N);
    
    fprintf('\n');
    
    fprintf('zero age lineages coalescence times:\n');
<<<<<<< HEAD
    fprintf('average: %f \n', mean_zero);
    fprintf('median: %f \n', median_zero);
=======
    fprintf('average: %g \n', mean_zero);
    fprintf('median: %g \n', median_zero);
    fprintf('number of simulations that reached coalescence: %g of %g \n', num_non_zero_sims_zero, iterations);
>>>>>>> L2Lv5
    
    fprintf('\n');
   
    fprintf('random age lineages coalescence times:\n');
<<<<<<< HEAD
    fprintf('average: %f \n', mean_random);
    fprintf('median: %f \n', median_random);    
    



=======
    fprintf('average: %g \n', mean_random);
    fprintf('median: %g \n', median_random);
    fprintf('number of simulations that reached coalescence: %g of %g \n', num_non_zero_sims_random, iterations);
>>>>>>> L2Lv5

fprintf('----------------------------------------------------\n');


<<<<<<< HEAD
function [number_generations, age_dist_m] = adjust_age_dist_m(burn_cycle_age_dist_m, age_dist_m)
    

    sizeof_age_dist_m = size(age_dist_m, 2);
    sizeof_burn_cycle_age_dist_m = size(burn_cycle_age_dist_m, 2);

    
    % did we make it out of the burn cycle?
    for i = 1:sizeof_burn_cycle_age_dist_m

        current_slice = burn_cycle_age_dist_m(:, i);
        next_slice = burn_cycle_age_dist_m(:, i + 1);
        
        if (isequal(current_slice, next_slice)) 
            number_generations = i;
            age_dist_m = burn_cycle_age_dist_m(:,1:i);
            fprintf("burn cycle too long: age_dist_m truncated at t = %d due to population stabilization.\n", i)
            return;
        end
    end
    
    for i = 1:sizeof_age_dist_m
        
        current_slice = age_dist_m(:, i);
        next_slice = age_dist_m(:, i + 1);
        
        if (isequal(current_slice, next_slice)) 
            number_generations = i;
            age_dist_m = age_dist_m(:,1:i);
            fprintf("age_dist_m truncated at t = %d due to population stabilization.\n", i)
            return;
        end
    end    
    % no adjustement required
    number_generations = sizeof_age_dist_m;
    age_dist_m = age_dist_m;
end
=======
>>>>>>> L2Lv5

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
    
<<<<<<< HEAD
end
=======
end

>>>>>>> L2Lv5
