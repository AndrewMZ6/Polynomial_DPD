function output_signal = generate_normalized_ofdm()

  # OFDM symbol check
  pkg load communications;

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


  fs = 50e6;
  fc = 10e6;
  t = 0:1/fs:(interpolated_size - 1)/fs;

  freqline = 0:fs/interpolated_size:fs - 1;

  Q_carr = -sin(2*pi*fc*t);
  I_carr = cos(2*pi*fc*t);

  rf_sig_ofdm = real(sig_ofdm_shifted_time).*I_carr + imag(sig_ofdm_shifted_time).*Q_carr;
  m = max(abs(rf_sig_ofdm));
  kkk = 1/m;
  rf_sig_ofdm = rf_sig_ofdm*kkk;
  output_signal = rf_sig_ofdm;

endfunction
