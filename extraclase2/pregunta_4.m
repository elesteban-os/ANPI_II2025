

function pregunta_4();
  clc; clear all; close all

  % Construir A y b
  A = zeros(10, 10);
  b = zeros(10, 1);

  for i = 1:10
    l = i - 1;
    sum = -1;
    for j = 1:10
      if (i == j)
        A(i, j) = i * 10;
        b(i, 1) = i;
        sum = 1;
        l = 1;
      else
        A(i, j) = l;
        l += sum;
      end
    endfor
  endfor

  % Aplicar eliminacion gaussiana
  display("Vector solucion: ");
  x = elimi_gauss(A, b)

endfunction


% ---- FUNC. MAIN ---
% x=elimi gauss(A,b), que resuelve un sistema de ecuaciones Ax = b, utilizando las funciones
% sust atras y triang sup
function x = elimi_gauss(A, b);
  % Calcular triangulo superior
  [At, bt] = triang_sup(A, b);

  % Obtener solucion de la ecuacion
  x = sust_atras(At, bt);

endfunction



% ---- FUNC. AUX ----
% [At,bt]=triang sup(A,b), que transforma un sistema Ax = b a un sistema Atx = bt usando
% operaciones elementales entre filas, y donde At es triangular superior.
function [At, bt] = triang_sup(A, b);
  % Obtener valor de m (cantidad de filas)
  m = size(A, 1);

  % Realizar matriz aumentada (A techo)
  At = [A b];

  % Algoritmo
  for k = 1:m-1
    for i = k + 1:m
      m_ik = At(i, k) / At (k, k);
      for j = k:m + 1
        At(i, j) = At(i, j) - m_ik * At(k, j);
      endfor
    endfor
  endfor

  bt = At(:, m + 1);
  At = At(:, 1:m);

endfunction


% ---- FUNC. AUX ----
% x=sust atras(A,b), que resuelve un sistema de ecuaciones Ax = b, donde A es triangular
% inferior, usando el metodo de sustitucion hacia atras.
function x = sust_atras(A, b);
  % Obtener valor de m (cantidad de filas)
  m = size(A, 1);
  x = zeros(1, m);

  % Recorrer las filas de las matrices de abajo para arriba
  for i = m:-1:1
    % Sumatoria de elementos anteriores
    sum = 0;
    for j = i + 1:m
      sum += A(i, j) * x(j);
    endfor

    % Calculo de una de las soluciones
    x(i) = (1 / A(i, i)) * (b(i) - sum);

  endfor

endfunction

