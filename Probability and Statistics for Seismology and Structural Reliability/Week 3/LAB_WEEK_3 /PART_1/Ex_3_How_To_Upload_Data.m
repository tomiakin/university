%% Initialisation of the workspace
clear % the "clear" command clear the workspace (imporant to clean the memory)
clc   % the "clc" clear the window command from previous command (only visual)


%% To upload data from a text file:
E = load('Data.txt');

% Data.txt is the name of the file containing the data
% E is the name of the vecotr/matrix contained in the txt file
% The symbol "=" assigne to "E" what is contained in the file 
% The command "load" does the job of taking the data from the file
% In the file only numbers are allowed;
% The ";" at the end avoids the vector to br printed in the "Command Window"

%% To upload data from an Excel file:
F = readtable('Data.xlsx',...
    'Sheet','Sheet2',...
    'Range','D3:E103');
G = table2array(F);

% Data.xlsx is the name of the file containing the data
% F is the name of the vecotr/matrix contained in the excel file
% The symbol "=" assigne to "F" what is contained in the file 
% The command "readtable" does the job of taking the data from the file
% The command "table2array" converts the table in a vector/matrix;
% The ";" at the end avoids the vector to br printed in the "Command Window"