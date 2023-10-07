clear
clc

%% Load the accelerogram accelerogram.txt
sampling_interval = 0.005;       % seconds
acc = load('accelerogram.txt');  % the unit is cm/s2

acc = acc / 100 / 9.81; % conversion of the unit to g

duration_accelerogram = (length(acc)-1)*sampling_interval;

time_vector = 0:sampling_interval:duration_accelerogram;

plot(time_vector,acc,'b')
xlabel('time [s]')
ylabel('acceleration [g]')
xlim([0 duration_accelerogram]);
ylim([-0.2 0.2])
grid on

%% Use the function response to calculate the response spectrum
vibration_period_range = 0:0.02:4;
damping = 0.05;
resp_spectrum = response(acc,sampling_interval,vibration_period_range,damping);

figure
plot(vibration_period_range,resp_spectrum,'b','linewidth',2)
xlabel('T - vibration period [s]')
ylabel('S_a - Spectral acceleration [g]')
xlim([0 4])
ylim([0 1.1*max(resp_spectrum)])
grid on


