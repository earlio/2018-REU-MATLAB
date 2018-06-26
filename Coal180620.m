
%% User inputs %%

% Po = input('Input an initial population vector and press enter:');
% numgen = input('Input a number of generations and press enter:');
% LM = input('Input a Leslie Matrix and press enter:');

%% Sample Inputs %% 

%SAMPLE 1
% Po = [2, 2, 2]; %sample initial population vector
% numgen = 10; %sample number of generations
% LM = [0 1 1.1; 0.6 0 0; 0 0.5 0]; %sample leslie matrix



%SAMPLE 2
Po = [2, 2, 2, 5]; %sample initial population vector
numgen = 20; %sample number of generations
LM = [0 1 1.1 1.2; 0.6 0 0 0; 0 0.5 0 0; 0 0 0.25 0]; %sample leslie matrix
life_table = [0 0.6 0; 1 0.5 1; 2 0.25 1.1; 3 0 1.12]

%SAMPLE 3
% Po = [56 72 16 84 23];
% numgen = 100; 
% LM = [0 0.5 0.75 1 1.25; 0.9 0 0 0 0; 0 0.75 0 0 0; 0 0 0.6 0 0; 0 0 0 0.4 0];

%SAMPLE 4
% Po = [506 823 234 348 294 1844]
% numgen = 50;
% LM = [0 0.1 0.1 0.1 0.1 0; 0.95 0 0 0 0 0; 0 0.92 0 0 0 0; 0 0 0.87 0 0 0; 0 0 0 0.56 0 0; 0 0 0 0 0.3 0];


%% Population Structure %%

t = 1:numgen; %create time vector from number of generations input
PopM = zeros(length(Po), length(t)); %creates a matrix which keeps track of the age structure of the population
PopM(:,1) = Po; %set the initial population structure equal to the given initial population
TotalPop = sum(Po); %determine the total population
Totals = zeros(1, length(t)); Totals(1) = TotalPop; %just as a test, see if the total population is maintained through time

for i = 2:length(t)
    PopM(:,i) = round(LM*PopM(:,i-1)); %Apply the leslie matrix to the previous population structure for each time step
    Totals(i) = sum(PopM(:,i)); %sums the total population of each generation (just for analysis)
end

cols = max(Totals); %sets the column dimension of the age structure matrix by determining the maximum population value derived from the Totals vector

AgeM = -1*ones(numgen,cols); %create a matrix with numgen rows and cols columns where every entry is negative one

for i = 1:length(t)
    AgeM(i,1:PopM(1,i)) = 0; %Boundary case, set the correct number of individuals to age 0
    for j = 2:length(Po)
        index = 1+sum(PopM(1:j-1,i)); %determine the number of individuals already assigned to an age, "+1" because MATLAB indices are inclusive
        AgeM(i,index:(index-1+PopM(j,i))) = j-1; %assign individuals of the next age group to the newest open spaces in the AgeM matrix
    end
end

AgeM %just to print AgeM
    
terminal_population = extract_terminal_population(AgeM);
 
disp(terminal_population);
% ez we'll start with just two hard-coded pairs to find the mrca of.

lineage_a_current_age = 0;
lineage_b_current_age = 0;


mrca = calculate_mrca(lineage_a_current_age, lineage_b_current_age, ...
    terminal_population, life_table, numgen);

disp(mrca)
