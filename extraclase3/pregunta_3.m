


function pregunta_3()
  pkg load symbolic
  clc; clear; close all

  % Parametros para llamar a la funcion
  yp = @(x, y) (x - y) / x;
  a = 2;
  b = 10;
  y0 = 4;
  n = 6;

  % Llamada a la funcion euler
  [x1, y1] = euler(yp, a, b, y0, n);

  % Llamada a la funcion predictor_corrector
  [x2, y2] = predictor_corrector(yp, a, b, y0, n);

  % Llamada a la funcion runge_kutta2
  [x3, y3] = runge_kutta2(yp, a, b, y0, n);

  % Llamada a la funcion runge_kutta3
  [x4, y4] = runge_kutta3(yp, a, b, y0, n);

  % Llamada a la funcion runge_kutta4
  [x5, y5] = runge_kutta4(yp, a, b, y0, n);

  % Llamada a la funcion Taylor 2
  syms xs ys;
  % usar la misma ecuacion simbolica que la version numerica: (x - y)/x
  yp_sym = (xs - ys) / xs;
  [x6, y6] = taylor2(yp_sym, a, b, y0, n);

  % Llamada a la funcion Adams Bashforth 2
  [x7, y7] = adams_bashforth2(yp, a, b, y0, n);

  % Llamada a la funcion Adams Bashforth 3
  [x8, y8] = adams_bashforth3(yp, yp_sym, a, b, y0, n);

  % Llamada a la funcion Adams Bashforth 4
  [x9, y9] = adams_bashforth4(yp, yp_sym, a, b, y0, n);


  % --- Gráfica ---
  figure
  hold on; grid on;
  plot(x1, y1, '-o', 'LineWidth', 1.8, 'MarkerSize', 6, 'MarkerFaceColor', 'b', 'DisplayName', 'Euler');
  plot(x2, y2, '-o', 'LineWidth', 1.8, 'MarkerSize', 6, 'MarkerFaceColor', 'r', 'DisplayName', 'Predictor-Corrector');
  plot(x3, y3, '-o', 'LineWidth', 1.8, 'MarkerSize', 6, 'MarkerFaceColor', 'y', 'DisplayName', 'Runge Kutta 2');
  plot(x4, y4, '-o', 'LineWidth', 1.8, 'MarkerSize', 6, 'MarkerFaceColor', 'g', 'DisplayName', 'Runge Kutta 3');
  plot(x5, y5, '-o', 'LineWidth', 1.8, 'MarkerSize', 6, 'MarkerFaceColor', "magenta", 'DisplayName', 'Runge Kutta 4');
  plot(x6, y6, '-o', 'LineWidth', 1.8, 'MarkerSize', 6, 'MarkerFaceColor', "cyan", 'DisplayName', 'Taylor 2');
  plot(x7, y7, '-o', 'LineWidth', 1.8, 'MarkerSize', 6, 'MarkerFaceColor', [0.5,0.2,0.8], 'DisplayName', 'Adams Bashforth 2');
  plot(x8, y8, '-o', 'LineWidth', 1.8, 'MarkerSize', 6, 'MarkerFaceColor', [0.2,0.6,0.4], 'DisplayName', 'Adams Bashforth 3');
  plot(x9, y9, '-o', 'LineWidth', 1.8, 'MarkerSize', 6, 'MarkerFaceColor', [0.9,0.4,0.1], 'DisplayName', 'Adams Bashforth 4');

  % Personalización del gráfico
  xlabel('x', 'FontSize', 12, 'FontWeight', 'bold');
  ylabel('y', 'FontSize', 12, 'FontWeight', 'bold');
  title('Comparación de métodos ODE', 'FontSize', 14, 'FontWeight', 'bold');
  legend('Location', 'northwest');
  axis tight;
  hold off;

endfunction

function [x, y] = euler(f, a, b, y0, n)
  h = (b - a) / n;
  x = a:h:b;
  y = zeros(1, n+1);
  y(1) = y0;
  for k=1:n
    y(k + 1) = y(k) + h * f(x(k), y(k));
  endfor
endfunction

function [x, y] = predictor_corrector(f, a, b, y0, n)
  h = (b - a) / n;
  x = a:h:b;
  y = zeros(1, n+1);
  y(1) = y0;

  for k=1:n
    K1 = f(x(k), y(k));
    K2 = f(x(k) + h, y(k) + h * K1);
    y(k + 1) = y(k) + h * (K1 + K2) * 0.5;
  endfor
endfunction

function [x, y] = runge_kutta2(f, a, b, y0, n)
  h = (b - a) / n;
  x = a:h:b;
  y = zeros(1, n+1);
  y(1) = y0;

  for k=1:n
    K1 = f(x(k), y(k));
    K2 = f(x(k) + (h / 2), y(k) + (h / 2) * K1);
    y(k + 1) = y(k) + h * K2;
  endfor
endfunction

function [x, y] = runge_kutta3(f, a, b, y0, n)
  h = (b - a) / n;
  x = a:h:b;
  y = zeros(1, n+1);
  y(1) = y0;

  for k=1:n
    K1 = f(x(k), y(k));
    K2 = f(x(k) + (h / 2), y(k) + (h / 2) * K1);
    K3 = f(x(k) + h, y(k) - h * K1 + 2 * h * K2);
    y(k + 1) = y(k) + (h / 6) * (K1 + 4 * K2 + K3);
  endfor
endfunction

function [x, y] = runge_kutta4(f, a, b, y0, n)
  h = (b - a) / n;
  x = a:h:b;
  y = zeros(1, n+1);
  y(1) = y0;

  for k=1:n
    K1 = f(x(k), y(k));
    K2 = f(x(k) + (h / 2), y(k) + (h / 2) * K1);
    K3 = f(x(k) + (h / 2), y(k) + (h / 2) * K2);
    K4 = f(x(k) + h, y(k) + h * K3);
    % falta
    y(k + 1) = y(k) + (h / 6) * (K1 + 2 * K2 + 2 * K3 + K4);
  endfor
endfunction

function [x, y] = taylor2(f, a, b, y0, n)
  h = (b - a) / n;
  x = a:h:b;
  y = zeros(1, n+1);
  y(1) = y0;

  syms xs ys;
  % f se espera como expresion simbolica en variables (xs, ys)
  dfdx = diff(f, xs);        % derivada parcial respecto a xs
  dfdy = diff(f, ys);        % derivada parcial respecto a ys
  fp_sym = simplify(dfdx + dfdy * f);

  % Convertir a funciones numericas
  f_num = matlabFunction(f, 'Vars', [xs, ys]);
  fp_num = matlabFunction(fp_sym, 'Vars', [xs, ys]);

  for k=1:n
    y(k + 1) = y(k) + h * f_num(x(k), y(k)) + (h^2 / 2) * fp_num(x(k), y(k));
  endfor
endfunction

function [x, y] = adams_bashforth2(f, a, b, y0, n)
  h = (b - a) / n;
  x = a:h:b;
  y = zeros(1, n+1);
  y(1) = y0;

  % Obtener paso inicial con Euler
  y(2) = y(1) + h * f(x(1), y(1));

  % Iterar
  for k=2:n
    y(k+1) = y(k) + (h/2) * (3*f(x(k), y(k)) - f(x(k-1), y(k-1)));
  endfor

endfunction

function [x, y] = adams_bashforth3(f, f_sym, a, b, y0, n)
  h = (b - a) / n;
  x = a:h:b;
  y = zeros(1, n+1);
  y(1) = y0;

  % Obtener dos pasos iniciales con Taylor
  [~, y_taylor] = taylor2(f_sym, a, a+2*h, y0, 2);
  y(2:3) = y_taylor(2:3);

  % Iterar
  for k=3:n
    y(k+1) = y(k) + (h/12) * (23*f(x(k), y(k)) - 16*f(x(k-1), y(k-1)) + 5*f(x(k-2), y(k-2)));
  endfor

endfunction

function [x, y] = adams_bashforth4(f, f_sym, a, b, y0, n)
  h = (b - a) / n;
  x = a:h:b;
  y = zeros(1, n+1);
  y(1) = y0;

  % Obtener tres pasos iniciales con Taylor
  [~, y_taylor] = taylor2(f_sym, a, a+3*h, y0, 3);
  y(2:4) = y_taylor(2:4);

  % Iterar
  for k=4:n
    y(k+1) = y(k) + (h/24) * (55*f(x(k), y(k)) - 59*f(x(k-1), y(k-1)) + 37*f(x(k-2), y(k-2)) - 9*f(x(k-3), y(k-3)));
  endfor

endfunction





