pkg load communications;
clearvars; close all; clc;

t = 0:0.001:pi;
f = 10;
sig = 0.6*sin(2*pi*f*t);

# sig2 = qammod(randint(1, 1000, [0, 15]), 16);


# OFDM

pkg load communications;
clearvars; close all; clc;

guards = 100;
fft_size = 1024;
interpolated_size = fft_size*20;
sig_ofdm = zeros(1, fft_size);
sc_num = fft_size - guards*2;
sig2 = qammod(randint(1, sc_num, [0, 15]), 16);

sig_ofdm(guards + 1:fft_size/2) = sig2(1:sc_num/2);
sig_ofdm(fft_size/2 + 1) = complex(0, 0);
sig_ofdm(fft_size/2 + 2:end-guards + 1) = sig2(sc_num/2 + 1:end);

sig_ofdm_shifted = fftshift(sig_ofdm);
sig_ofdm_shifted = [sig_ofdm_shifted(1:fft_size/2), complex(zeros(1, interpolated_size - fft_size)), sig_ofdm_shifted(fft_size/2 + 1: end)];

sig_ofdm_shifted_time = ifft(sig_ofdm_shifted);


figure;
plot(abs(fftshift(sig_ofdm)));

scatterplot(sig_ofdm);


fs = 100e6;
fc = 10e6;
t = 0:1/fs:(interpolated_size - 1)/fs;

freqline = 0:fs/interpolated_size:fs - 1;

% TX

Q_carr = -sin(2*pi*fc*t);
I_carr = cos(2*pi*fc*t);

rf_sig_ofdm = real(sig_ofdm_shifted_time).*I_carr + imag(sig_ofdm_shifted_time).*Q_carr;
m = max(abs(rf_sig_ofdm));
k = 1/m;
rf_sig_ofdm = rf_sig_ofdm*k;

figure;
plot(freqline, abs(fft(rf_sig_ofdm)));



% PA model
A = 1.2;
B = 1.43;
r = -1:0.1:1;

PA_saleh = @(r) (A*r)./(1 + B*(r.^2));


% PA model end

rf_sig_ofdm = PA_saleh(rf_sig_ofdm);


figure;
plot(freqline, abs(fft(rf_sig_ofdm)));

% RX

scatterplot(fft(rf_sig_ofdm));



y = PA_saleh(r);


figure;
plot(r, y);
grid on;

return;

figure;
plot(t, sig);
grid on;

figure;
plot(t, PA_saleh(sig));
grid on;



scatterplot(PA_saleh(sig2));


max(abs(sig2))
