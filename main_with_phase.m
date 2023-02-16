close all; clearvars; clc;

A_array = 0.5:0.3:2;
B_array = 0.5:0.3:2;
A_array(end) = 1.2;
B_array(end) = 1;







% PA model
A = A_array(end);
B = B_array(end);
alpha = 0.37;
beta = 0.68;

PA_saleh = @(r) (A*r)./(1 + B*(r.^2));
PA_saleh_phase = @(r) (alpha*(r.^2))./(1 + beta*(r.^2));

#a1 = 0.3; a2 = 0.6; a3 = 0.2;
#PA_saleh_phase = @(r) a1*(r.^2) + a2*(r.^3) - a3*(r.^5);


x = linspace(-1, 1, 100);

figure;
plot(x, PA_saleh_phase(x));


#return;

am_array = PA_saleh(x);
pm_array = PA_saleh_phase(x);
saleh_result = am_array.*(exp(i*pm_array));


PA_model_order = 5;

for i = 1:PA_model_order + 1
  H(:, i) = x'.^(i-1);
end

y = PA_saleh(x)';
maxy = max(abs(y));
y = y./maxy;

PA_model_Theta = inv(H'*H)*H'*y;

model = zeros(1, length(x));

figure;
for j = 1:PA_model_order + 1
  # На каждой итерации возрастает степень полинома
  # Все выходные значения для текущего полинома вычисляются одновременно
  mm = plot(x, y, 'LineWidth', 2); hold on;
  ylim([-1, 1]);
  model = model + PA_model_Theta(j)*x.^(j-1);
  h = plot(x, model, 'LineWidth', 1, 'r+');
  grid on;

  legend('PA', 'model');
  pause(0.05);
  clear h;
  hold off;
end



newx = model;

DPD_poly_order = 11;
newy = x;

for i = 1:DPD_poly_order + 1
  newH(:, i) = newx'.^(i-1);
end



DPD_Theta = inv(newH'*newH)*newH'*newy';

figure;

invmodel = zeros(1, length(newx));

for j = 1:DPD_poly_order + 1
  mm = plot(newx, newy, 'LineWidth', 2); hold on;
  ylim([-1, 1]);
  xlim([-1, 1]);
  invmodel = invmodel + DPD_Theta(j)*newx.^(j-1);
  h = plot(newx, invmodel, 'LineWidth', 1, 'r+');
  grid on;

  legend('inverse PA', 'DPD model');
  pause(0.05);
  clear h;
  hold off;
end


#PA_model = @(x) Theta(1) + Theta(2)*x.^1 + Theta(3)*x.^2 + Theta(4)*x.^3 + Theta(5)*x.^4 + Theta(6)*x.^5 + Theta(7)*x.^6 + Theta(8)*x.^7;
#DPD_model = @(x) newTheta(1) + newTheta(2)*x.^1 + newTheta(3)*x.^2 + newTheta(4)*x.^3 + newTheta(5)*x.^4 + newTheta(6)*x.^5 + newTheta(7)*x.^6 + newTheta(8)*x.^7 + ...
#                   newTheta(9)*x.^8 + newTheta(10)*x.^9 + newTheta(11)*x.^10 + newTheta(12)*x.^11;


# Синусоида для тестирования обратного полинома

fc1 = 1e3;  A = 1; fs = 10e4; N = 10000; fc2 = 2e3;
tsin = 0:1/fs:(N-1)/fs;
sig = sin(2*pi*fc1*tsin) + sin(2*pi*fc2*tsin) ;
sigmax = max(sig);
sig = sig./sigmax;
sig = sig*A;


after_PA_model = _apply_polynomial(PA_model_Theta, sig);
after_PA_model_with_DPD = _apply_polynomial(PA_model_Theta, _apply_polynomial(DPD_Theta, sig));


figure;
plot(tsin, after_PA_model, tsin, after_PA_model_with_DPD);
grid on;
legend('without DPD', 'with DPD');


figure;
plot(abs(fft(after_PA_model)), 'LineWidth', 2); hold on;
plot(abs(fft(after_PA_model_with_DPD)), 'LineWidth', 2.5); hold off;
legend('without DPD', 'with DPD');
grid on;




figure;

plot(x, _apply_polynomial(PA_model_Theta, x), 'LineWidth', 2); hold on;
plot(newx, _apply_polynomial(DPD_Theta, newx), 'LineWidth', 2); hold on;
plot(x, _apply_polynomial(PA_model_Theta, _apply_polynomial(DPD_Theta, x)), 'LineWidth', 2); hold off;
xlim([-1, 1]);
ylim([-1, 1]);
grid on;
legend('PA AM-AM', 'DPD AM-AM', 'PA + DPD AM-AM');



figure;
plot(x, PA_saleh(_apply_polynomial(DPD_Theta,x)));
grid on;
title('DPD + PA AM-AM');




pkg load communications;
#clearvars; close all; clc;




guards = 100;
fft_size = 1024;
interpolated_size = fft_size*20;
sig_ofdm = zeros(1, fft_size);
sc_num = fft_size - guards*2;
sig2 = qammod(randint(1, sc_num, [0, 63]), 64);

sig_ofdm(guards + 1:fft_size/2) = sig2(1:sc_num/2);
sig_ofdm(fft_size/2 + 1) = complex(0, 0);
sig_ofdm(fft_size/2 + 2:end-guards + 1) = sig2(sc_num/2 + 1:end);

sig_ofdm_shifted = fftshift(sig_ofdm);
#sig_ofdm_shifted = sig_ofdm_shifted.*exp(i*pm_array);
sig_ofdm_shifted = [sig_ofdm_shifted(1:fft_size/2), complex(zeros(1, interpolated_size - fft_size)), sig_ofdm_shifted(fft_size/2 + 1: end)];

sig_ofdm_shifted_time = ifft(sig_ofdm_shifted);

scatterplot(fft(sig_ofdm_shifted_time))
#####
i = sqrt(-1);


figure;
for p = 0.1:0.01:0.99

  gg = abs(sig_ofdm_shifted);
  angles = angle(sig_ofdm_shifted);
  gg = (gg./max(gg))*p;
  max(gg)
  am_array = PA_saleh(_apply_polynomial(DPD_Theta,gg));
  pm_array = PA_saleh_phase(gg);
  max(pm_array)
  saleh_result = am_array.*(exp(i*(pm_array + angles))); #
  ##
  ##

  kk = scatter(real(saleh_result), imag(saleh_result), 'filled');
  grid on;
  xlim([-1, 1]);
  ylim([-1, 1]);
  pause(0.01);
  clear kk;
end


figure;
plot((pm_array + angles)./angles);

return;
###


fs = 50e6;
fc = 10e6;
t = 0:1/fs:(interpolated_size - 1)/fs;

freqline = 0:fs/interpolated_size:fs - 1;

% TX

Q_carr = -sin(2*pi*fc*t);
I_carr = cos(2*pi*fc*t);

rf_sig_ofdm = real(sig_ofdm_shifted_time).*I_carr + imag(sig_ofdm_shifted_time).*Q_carr;
m = max(abs(rf_sig_ofdm));
k = 1/m;
rf_sig_ofdm = rf_sig_ofdm*k*A;
max(abs(rf_sig_ofdm))


figure;
plot(rf_sig_ofdm);


after_PA_model = _apply_polynomial(PA_model_Theta, rf_sig_ofdm);
after_PA_model_with_DPD = _apply_polynomial(PA_model_Theta, _apply_polynomial(DPD_Theta, rf_sig_ofdm));

figure;
plot(freqline, 10*log10(abs(fft(after_PA_model)))); hold on;
plot(freqline, 10*log10(abs(fft(after_PA_model_with_DPD))));
title('DPD + polynomial');
grid on;

scatterplot(fft(after_PA_model));
scatterplot(fft(after_PA_model_with_DPD));


after_saleh_model = PA_saleh(rf_sig_ofdm);
after_saleh_model_with_DPD = PA_saleh(_apply_polynomial(DPD_Theta, rf_sig_ofdm));


figure;
plot(freqline, 10*log10(abs(fft(after_saleh_model)))); hold on;
plot(freqline, 10*log10(abs(fft(after_saleh_model_with_DPD))));
title('DPD + saleh');
grid on;



scatterplot(fft(after_saleh_model));
scatterplot(fft(after_saleh_model_with_DPD));
