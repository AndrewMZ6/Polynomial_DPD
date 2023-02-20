close all; clearvars; clc;
pkg load communications;


% PA model



k = 1;
sc_num = 1000;
sig = qammod(randint(1, sc_num, [0, 63]), 64);
normalized_sig = sig/max(abs(sig))*k;

amplified_sig = saleh_PA(scaled_sig);



