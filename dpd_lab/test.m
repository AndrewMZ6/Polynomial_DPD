clearvars; close all; clc;
pkg load communications;


l = 0.7;
x = -l:0.001:l;
y = amplifier_model(x);


[y1, x1] = max(y)
uu = x(x1)

figure;
plot(x, y);
grid on;
