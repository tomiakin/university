%% Initialisation of the workspace
clear % the "clear" command clear the workspace (imporant to clean the memory)
clc   % the "clc" clear the window command from previous command (only visual)

load SHAREver3

ID = 1:1:length(SHAREver3.Mw);

Year = SHAREver3.Year;
Month = SHAREver3.Mo;
Day = SHAREver3.Da;
Hour = SHAREver3.Ho;
Minute = SHAREver3.Mi;
Second = SHAREver3.Se;

Earthquake_date = datetime(Year,Month,Day,Hour,Minute,Second);

plot(sort(Earthquake_date),ID)
xlabel('Time')
ylabel('Even number')