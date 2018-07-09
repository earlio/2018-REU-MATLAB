function [leslie_matrix,ages,scaling, life_table_m] = life_to_leslie(file_path_name,rescale)

% function [leslie_matrix,ages,scaling] = life_to_leslie()

% A function to input a life table and produce a Leslie matrix that can be used to simulate population demography over time. 

% Inputs:
% a life table in a .csv or an .xlsx file in the form of an number of ages x 3 matrix. 
% Column 1 contains the ages as integers, Column 2 contains age-specific survival as a real number, Column 3 contains age-specific fecundity as a real number. 
% Life tables should include a final age class where both survival and fecundity are zero.
%
% file_path_name - path name to life table file
% rescale - boolean for whether to rescale Leslie matrix such that populaiton size is constant over time

% Outputs:
%
% leslie_matrix - corresponding Leslie Matrix that is square ages x ages
% ages - integer number of age classes in life table
% scaling - real number pop growth rate of original life table


    %filename = file; %assign the file to a variable
    
    %file_path_name = '/Users/matthewhamilton/Downloads/reu2018gu-es-p-master/Sample_LT1.xlsx';
    
    % use importdata to read life table file since it does not require
    % knowledge of file type, works for Excel and for .csv
    file_contents = importdata(file_path_name);
    
    life_table_m = file_contents.data;
    
%     % read in life tabel file if we know file is either csv or .xls
%     if isequal(filename(end-2:end),'csv')
%         life_table_m = csvread(filename,2,0); %read the csv file and assign it to a matrix called life_table_m. 2,0 offsets the entries to begin with teh values
%     elseif isequal(filename(end-2:end),'lsx')
%         life_table_m = xlsread(filename, cell_range);
%     end

    [ages] = size(life_table_m,1); %determine the number of rows of the life_table to create a correctly sized leslie matrix

    leslie_matrix = zeros(ages,ages); %create a square leslie matrix with dimension = number of ages in the life table

    for i = 1:ages-1 
       
        leslie_matrix(1,i) = life_table_m(i,3); %set the first row of the leslie matrix to the age-specific fecundity coefficients
       
        leslie_matrix(i+1,i) = life_table_m(i,2); %set the diagonal (beginning at (2,1)) equal to the age-specific survival coefficients
    end

    % get eigenvalues of Leslie matrix, leading eigenvalue is equal to pop growth rate
    lambda = eig(leslie_matrix); 

    if lambda(1) > 0
        scaling = lambda(1);
    else
        scaling = lambda(2);
    end

    % test of rescaling
    %rescale = true;
    
    if rescale == true
        leslie_matrix = leslie_matrix./scaling;
    end


