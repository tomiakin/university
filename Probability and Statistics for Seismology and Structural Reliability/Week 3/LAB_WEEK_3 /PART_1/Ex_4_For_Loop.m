%% Initialisation of the workspace
clear % the "clear" command clear the workspace (imporant to clean the memory)
clc   % the "clc" clear the window command from previous command (only visual)


%% Definition of a FOR loop
% The for loop is a command that allows to repete a command or a list of
% commands in a cyclic manner.

N = 100;
for i = 1:N
    
    % Here you write the instruction you want to run
    disp(num2str(i)) % this is the easier command possible
    
end

% In the example provided above "N" is the number of times that the cycle
% is performed. Specifically, from 1 to 100 the cycle will print on the
% command window the progressiv vale "i".

% Obviously, this is the easiest command possible. In the for loop you can
% have several command lines.

%% Definition of a WHILE loop
% The while loop is a command that allows to repete a command or a list of
% commands in a cyclic manner until a certain condition is obtained.

N = 100;
i = 1;
while i<=N
    
    % Here you write the instruction you want to run
    disp(num2str(i)) % this is the easier command possible
    i = i+1;
    
end

% In the example above, I created a while loop that does the same
% instructions of the for loop presented earlier.