close all; clearvars; clc;


% PA model
A = 1.2;
B = 1;
alpha = 0.37;
beta = 0.68;

PA_saleh = @(r) ((A*r)./(1 + B*(r.^2)))./0.6;
PA_saleh_phase = @(r) ((alpha*(r.^2))./(1 + beta*(r.^2)));


x = linspace(-1, 1, 100);
x_abs = abs(x);


am_array = PA_saleh(x);
pm_array = PA_saleh_phase(x_abs);
pm_max = max(pm_array);


k = 1;
sc_num = 1000;
sig = qammod(randint(1, sc_num, [0, 15]), 16);
sig = sig/max(abs(sig))*k;

amps = abs(sig);
angles = angle(sig);


figure;

upperlim = 1;
z = zeros(1, length(x));
z2 = zeros(1, length(x));
pause(3);
for p = 0.1:0.01:upperlim
  iter_sig = sig*p;

  amps = abs(iter_sig);
  angles = angle(iter_sig);
  new_amps = PA_saleh(amps);
  add_phases = PA_saleh_phase(amps);

  new_sig = new_amps.*exp(i*(angles + add_phases));
  new_sig2 = new_amps.*exp(i*(angles));


  subplot(2, 2, [3, 4]);
    h = scatter(real(new_sig), imag(new_sig), 'filled'); hold on;
    #l2 = scatter(real(iter_sig), imag(iter_sig), 'filled'); hold off;
    l = scatter(real(new_sig2), imag(new_sig2), 'filled'); hold off;

    lims = [-upperlim, upperlim];
    xlim(lims);
    ylim(lims);
    legend('amplifier with phase', 'amplifier without phase');
    grid on;

  subplot(2, 2, 2);
    z(:) = p;
    z2(:) = PA_saleh(p);
    o = plot(x, am_array); hold on;
    o2 = plot(z, x); hold on;
    o3 = plot(x, z2); hold on;
    o4 = plot(p, PA_saleh(p), 'r+', 'LineWidth', 4); hold off;
    grid on;
    xlabel('PA input');
    ylabel('PA output');


  subplot(2, 2, 1);
    z2(:) = PA_saleh_phase(p);
    x_abs2 = linspace(0, 0.25, length(z));
    o = plot(x_abs, pm_array); hold on;
    o2 = plot(z, x_abs2); hold on;
    o3 = plot(x_abs, z2); hold on;
    o4 = plot(p, PA_saleh_phase(p), 'r+', 'LineWidth', 4); hold off;
    grid on;
    xlabel('PA input');
    ylabel('PA phase');


  pause(0.01);
  clear h;
  clear o; clear o2; clear o3; clear o4;
  clear l;

end
