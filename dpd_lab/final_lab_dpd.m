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


dpd_values = zeros(1, length(input_values));

for i = 1:length(input_values)
    dpd_values(i) = M(input_values(i));
    disp(['Filling dpd_values. Iteration ', num2str(i)]);
endfor

plot(input_values, amplifier_model(dpd_values));


rf_sig_ofdm = generate_normalized_ofdm();
rf_sig_ofdm = rf_sig_ofdm*0.316;


# This is comment line from Polynomial DPD repo
# newGit COMMENT! 

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
  legend('No DPD', 'DPD');


a1 = fft(amplified_rf_sig_ofdm);
a2 = fft(amplifier_model(dpd_rf_sig_ofdm));

a1_cut = a1(4097-512:4097+512);
a2_cut = a2(4097-512:4097+512);


scatterplot(a1_cut);
title('No DPD');
scatterplot(a2_cut);
title('DPD');



