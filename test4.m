# comm
clearvars; close all; clc;


# input values
N = 100; r = 10;
x = rand(1, N)*2 - 1;
x = linspace(0, r, N);


K = 7;
M = 5;

oo = 100;


figure;
  while oo > 0
    y = zeros(1, N);
    a = rand(K, M + 1);#*randint(1, 1);

    for n = 1:N
      for k = 1:K
        for m = 0:M

          temp = n - m;

          if temp > 0
            y(n) = y(n) + a(k, m + 1)*x(temp)*(abs(x(temp)))^(k-1);
          endif

        endfor
      endfor
    endfor

    plot(x, y); hold on;
    grid on;

    pause(0.1);
    oo = oo - 1;
  endwhile



