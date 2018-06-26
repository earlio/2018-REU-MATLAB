function [terminal_population] = extract_terminal_population2(population_matrix)
%returns the portion of the last row which contains individuals
    [ages] = size(population_matrix,1); %determines the number of rows of the population matrix to extract the number of age cohorts
    terminal_population = []; %initializes the terminal population vector
    for i=1:ages
        terminal_population = [terminal_population, ((i-1)*ones(1,population_matrix(i,end)))]; %concatenates the terminal population vector with a vector corresponding to each age cohort
    end
end