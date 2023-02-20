close all; clearvars; clc;


x = linspace(0, 1, 100);
phases = linspace(0, 4*pi, 100);

test_sig = x.*exp(i*(phases));
test_sig_PA = saleh_PA(test_sig);
max(abs(test_sig))


before_amp_color = '#990000';
after_amp_color = '#004080';

figure;
  subplot(1, 2, 1);
    plot(real(test_sig), imag(test_sig), 'color', before_amp_color, 'marker', 'o', 'markersize', 5, 'LineStyle', 'None'); hold on;
    plot(real(test_sig_PA), imag(test_sig_PA), 'color', after_amp_color, 'marker', 'd', 'markersize', 5, 'LineStyle', 'None'); hold on;
    grid on;

    for _ = 1:length(test_sig)
      plot([real(test_sig(_)), real(test_sig_PA(_))], [imag(test_sig(_)), imag(test_sig_PA(_))], 'color', '#336699', 'LineStyle', '-'); hold on;

    endfor

  subplot(1, 2, 2);
    plot(x, abs(test_sig), 'color', before_amp_color, 'marker', 'o' , 'markersize', 5, 'LineStyle', '-'); hold on;
    plot(x, abs(test_sig_PA), 'color', after_amp_color, 'marker', 'd' , 'markersize', 5, 'LineStyle', '-'); hold off;
    grid on;



