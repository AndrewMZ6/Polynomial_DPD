function output_signal = saleh_PA(input_signal)
  A = 1.2;
  B = 1;
  alpha = 0.37;
  beta = 0.68;

  saleh_AM = @(r) ((A*r)./(1 + B*(r.^2)));
  saleh_PH = @(r) (alpha*(r.^2))./(1 + beta*(r.^2));

  initial_phases = angle(input_signal);
  initial_modules = abs(input_signal);
  m = max(initial_modules);

  if m > 1
    warning('saleh_PA function -> input signal module is > 1. Make sure to normalize it!');
  endif

  new_amplitudes = saleh_AM(initial_modules);
  additional_phases = saleh_PH(initial_modules);

  output_signal = new_amplitudes.*exp(i*(initial_phases + additional_phases));

endfunction
