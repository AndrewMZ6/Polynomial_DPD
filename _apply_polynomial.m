function result = _apply_polynomial(Theta, input_data)

  # Порядок полинома всегда на 1 меньше количества элементов в Theta
  K = length(Theta);

  result = zeros(1, length(input_data));

  for i = 1:K
    result = result + Theta(i)*input_data.^(i-1);
  endfor

end
