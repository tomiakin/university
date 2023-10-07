%% Initialisation of the workspace
clear % the "clear" command clear the workspace (imporant to clean the memory)
clc   % the "clc" clear the window command from previous command (only visual)


%% To create a Matrix C (3x6):
C = [1 2 3 4 5 6; 7 8 7 8 7 8; 1 2 3 4 5 6];

% C is the name of the vecotr
% The symbol "=" assigne to "C" what is defined on the right-hand-side 
% The "[]" parenthesis define the matrix
% Aving numbers in order in the parentesis witout commas define the raws;
% the semicolumns at the end of each columns tell matlab when to shift to
% the new line.
% The ";" at the end avoids the vector to br printed in the "Command Window"

%% To create a Matrix D (3x6) in a more clear manner:
D = [1 2 3 4 5 6; ...
     7 8 7 8 7 8; ...
     1 2 3 4 5 6];

% D is the name of the vecotr
% The symbol "=" assigne to "D" what is defined on the right-hand-side 
% The "[]" parenthesis define the vector
% Aving numbers in order in the parentesis witout commas define the raws
% The symbol "..." tells to matlab that the command continues to the matrix
% line.
% The ";" at the end avoids the vector to br printed in the "Command Window"