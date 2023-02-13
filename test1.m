close all; clearvars; clc;

A_array = 0.5:0.3:2;
B_array = 0.5:0.3:2;
A_array(end) = 1.2;
B_array(end) = 1.43;

figure;

for p = 1:length(A_array)
  % PA model
  #A = 1.2;
  A = A_array(p);
  #B = 1.43;
  B = B_array(p);
  #r = -1:0.1:1;

  PA_saleh = @(r) (A*r)./(1 + B*(r.^2));


  x = -1:0.05:1;

  K = 7;

  for i = 1:K + 1
    x_matrix(:, i) = x'.^(i-1);
  end

  y_matrix = PA_saleh(x)';

  coeffs = inv(x_matrix'*x_matrix)*x_matrix'*y_matrix;
  coeffs2 = inv(y_matrix'*y_matrix)*y_matrix'*x_matrix;

  #newf = @(x) coeffs(1) + coeffs(2)*x + coeffs(3)*x.^2 + coeffs(4)*x.^3 + coeffs(5)*x.^4 + coeffs(6)*x.^5;
  #if p == 1
  newf2 = zeros(1, length(x));
  newf3 = zeros(1, length(x));
  #end

  for j = 1:K + 1
    mm = plot(x, y_matrix, 'LineWidth', 2); hold on;
    ylim([-1, 1]);
    newf2 = newf2 + coeffs(j)*x.^(j-1);
    newf3 = newf3 + coeffs2(j)*y_matrix.^(j-1);
    h = plot(x, newf2, 'LineWidth', 1, 'r+');
    grid on;

    legend('PA', 'model');
    pause(0.05);
    clear h;
    hold off;
  end

pause(2);
clear mm;

end

#figure;
#plot(y_matrix, newf3);

newf = @(x) coeffs2(1) + coeffs2(2)*x + coeffs2(3)*x.^2 + coeffs2(4)*x.^3 + coeffs2(5)*x.^4 + coeffs2(6)*x.^5;


