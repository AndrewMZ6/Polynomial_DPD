close all; clearvars; clc;
pkg load communications;


% PA model



k = 1;
sc_num = 1000;
sig = qammod(randint(1, sc_num, [0, 63]), 64);
normalized_sig = sig/max(abs(sig))*k;

figure;

for p = 0.1:0.01:1

  scaled_sig = normalized_sig*p;
  amplified_sig = saleh_PA(scaled_sig);
  h = scatter(real(amplified_sig), imag(amplified_sig), 'filled');
  grid on;
  xlim([-1, 1]);
  ylim([-1, 1]);

  pause(0.1);
  clear h;

endfor

