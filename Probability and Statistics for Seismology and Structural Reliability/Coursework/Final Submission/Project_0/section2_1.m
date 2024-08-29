 close all
 clear all
 clc

u_m = 1843.43;
s_m = 150.16;
u_p = 1011.7;
s_p = 50.16;
P = 3055.14;
M = 4144.55;

pm = 0.91;

B1 = (1 - (u_m/M) - (u_p/P)) / (sqrt((s_m/M)^2 + (s_p/P)^2 + 2 * (1/M) * (1/P) * pm * s_m * s_p));

B2 = (M*P - P*u_m - M*u_p) / (sqrt((P*s_m)^2 + (M*s_p)^2 + 2 * P * M * pm * s_m * s_p));

