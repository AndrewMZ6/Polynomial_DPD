clearvars; close all; clc;
pkg load communications;

## Number of digits after decimal point
n = 3;
precision = 1/(10^n);


## Create input values array
input_values = -0.4:precision:0.4;
input_values = round2(input_values, n);

## Define desirable amplification
k = 11.489;


er = 0.00001;
d = 0.1;

# Number of while loop iterations until break
while_break_count = 1000;
for i = 1:length(input_values)
    e = 1; phi = 0;
    e_cache = 0;
    s = input_values(i);
    ppp = 1;
    while abs(e) > er
        if ppp < while_break_count
            A = s + phi;
            P = amplifier_model(A);
            P = P./k;
            e = s - P;
            if e_cache != e
                e_cache = e;
                phi = phi + d*e;
                ppp = ppp + 1;
            else
                disp(["e_cache break for ", num2str(i), "'th element after ", num2str(ppp), " iterations."]);
                break
                ppp = ppp + 1;
            endif
        else
            disp(['End of ', num2str(while_break_count) ,' iterations']);
            break
        endif

    end

    arr(i) = phi;
end


right_column = input_values + arr;

M = containers.Map("KeyType", "double", "ValueType", "double");
for i = 1:length(input_values)
    M(input_values(i)) = right_column(i);
    disp(['Filling hash table. Iteration ', num2str(i)]);
endfor

figure;
plot(input_values, amplifier_model(input_values)); grid on; hold on;


dpd_values = zeros(1, length(input_values));

for i = 1:length(input_values)
    dpd_values(i) = M(input_values(i));
    disp(['Filling dpd_values. Iteration ', num2str(i)]);
endfor

plot(input_values, amplifier_model(dpd_values)); hold off;

figure;
plot(input_values, amplifier_model(dpd_values));
grid on;


figure;
plot(input_values, dpd_values);




# OFDM symbol check
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


A = 0.316;

fs = 50e6;
fc = 10e6;
t = 0:1/fs:(interpolated_size - 1)/fs;

freqline = 0:fs/interpolated_size:fs - 1;

Q_carr = -sin(2*pi*fc*t);
I_carr = cos(2*pi*fc*t);

rf_sig_ofdm = real(sig_ofdm_shifted_time).*I_carr + imag(sig_ofdm_shifted_time).*Q_carr;
m = max(abs(rf_sig_ofdm));
kkk = 1/m;
rf_sig_ofdm = rf_sig_ofdm*kkk*A;


amplified_rf_sig_ofdm = amplifier_model(rf_sig_ofdm);

for i = 1:length(rf_sig_ofdm)
  temp = round2(rf_sig_ofdm(i), n);
  if (!temp)
    dpd_rf_sig_ofdm(i) = 0;
  else

    dpd_rf_sig_ofdm(i) = M(temp);
  endif
endfor


bla1 = abs(fft(amplified_rf_sig_ofdm));
m = min(bla1);
bla = abs(fft(amplifier_model(dpd_rf_sig_ofdm)));


zer = find(bla == 0);
bla(zer) = m;


figure;
  plot(10*log10(bla1)); hold on;
  plot(10*log10(bla)); hold off;


a1 = fft(amplified_rf_sig_ofdm);
a2 = fft(amplifier_model(dpd_rf_sig_ofdm));

a1_cut = a1(4097-512:4097+512);
a2_cut = a2(4097-512:4097+512);


scatterplot(a1_cut);
scatterplot(a2_cut);



