

function pregunta_2()
  clc; clear; close all

  # Valores de entrada
  f = @(x) log(asin(x)) ./ log(x);    # Por medio de propiedades de logaritmos, f(x) = log_x (sin^(−1)(x))
  a = 0.1;
  b = 0.8;
  xv = 0.1:0.1:0.8;
  x_0 = 0.55;

  # Calcular la cota de interpolación
  ct = cota_interpolacion(f, a, b, xv, x_0);
  disp(['Cota de interpolación en x_0 = ', num2str(x_0), ': ', num2str(ct)]);

end

function ct = cota_interpolacion(f, a, b, xv, x_0)
  pkg load symbolic
  syms x;

  # Obtener grado del polinomio
  n = length(xv) - 1;

  # Convertir a función simbólica
  f_sym = sym(f);

  # Calcular la n-ésima derivada
  f_n_deriv = diff(f_sym, x, n + 1);

  # Convertir de nuevo a función numérica
  f_n_deriv_func = matlabFunction(f_n_deriv);

  # Obtener el máximo valor absoluto de la n-ésima derivada en [a, b]
  f_n_deriv_max = max(abs(f_n_deriv_func(linspace(a, b, 2000))));

  # Calcular el producto (x_0 - x_0)(x_0 - x_1)...(x_0 - x_n)
  prod_term = 1;
  for i = 1:length(xv)
    prod_term = prod_term * (x_0 - xv(i));
  endfor

  # Calcular la cota de interpolación
  ct = (f_n_deriv_max / factorial(n + 1)) * abs(prod_term);
endfunction
