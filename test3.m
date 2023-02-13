close all; clearvars; clc;

A_array = 0.5:0.3:2;
B_array = 0.5:0.3:2;
A_array(end) = 1.2;
B_array(end) = 1.43;

figure;

% PA model
#A = 1.2;
A = A_array(end);
#B = 1.43;
B = B_array(end);
#r = -1:0.1:1;

PA_saleh = @(r) (A*r)./(1 + B*(r.^2));
#gg = [0.308429099596450 0.286349689341066 0.359712939604602 0.503475227374009 0.733917596608961 0.842448844089383 1.04830148276661 1.29587806964960];
gg = [0.0904507710557533 0.112083500602168 0.144254835039818 0.185413416536662 0.232080200501253 0.278392597968070 0.325256384476171 0.385256988277728 0.439142290530491 0.499520876112252 0.555854241338112 0.615471394037067 0.674263674614306 0.739957306423443 0.810816634346046];
maxgg = max(gg);
normgg = gg./maxgg;
length(gg)
xgg = linspace(0, 1, length(gg));

figure;
plot(xgg, normgg);

#x = -1:0.05:1;
x = xgg;

K = 5;

for i = 1:K + 1
  H(:, i) = x'.^(i-1);
end

#y = PA_saleh(x)';
y = normgg';

Theta = inv(H'*H)*H'*y;

newf2 = zeros(1, length(x));




for j = 1:K + 1
  # На каждой итерации возрастает степень полинома
  # Все выходные значения для текущего полинома вычисляются одновременно
  mm = plot(x, y, 'LineWidth', 2); hold on;
  ylim([-1, 1]);
  newf2 = newf2 + Theta(j)*x.^(j-1);
  h = plot(x, newf2, 'LineWidth', 1, 'r+');
  grid on;

  legend('PA', 'model');
  pause(0.1);
  clear h;
  hold off;
end

err = newf2 - normgg;
figure;
  plot(err);
  title('error');

std(err)

