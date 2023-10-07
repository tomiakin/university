%% Initialisation of the workspace
clear % the "clear" command clear the workspace (imporant to clean the memory)
clc   % the "clc" clear the window command from previous command (only visual)


%% To create a Structure A containing numbers, vector, and text:
A(1).Name = 'First';
A(1).X    = 3.14;
A(1).Y    = [1 2 3 4];
A(1).Z    = [1 2 3 4; ...
             5 6 7 8; ...
             9 8 7 6];

% A is the name of the structure
% A(1) indicates the first element in the structure
% A(1).Name = 'First' indicates that the filed Name of the structure A is the string 'First'
% A(1).X = 3.14 indicates that the element X of the structure A is the number 3.14
% A(1).Y = [1 2 3 4]; indicates that the element Y of the structur A is the vector [1 2 3 4]
% A(1).Z = ...; indicates that the element Y of the structur A is the matrix indicated above
% The ";" at the end avoids the vector to br printed in the "Command Window"

%% To add a new element to the Structure A:
A(2).Name = 'Home';
A(2).X    = 567;
A(2).Y    = [4 4 4 4];
A(2).Z    = [0 0 0 0; ...
             1 1 1 1; ...
             2 2 2 2];