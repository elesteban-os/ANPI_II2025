

function pregunta_4();
  ## Función principal que construye un sistema de ecuaciones Ax = b
  ## y lo resuelve usando eliminación gaussiana

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
function x = elimi_gauss(A, b);
  ## Resuelve un sistema de ecuaciones lineales Ax = b utilizando
  ## el método de eliminación gaussiana
  ##
  ## Parámetros:
  ##   A    matriz de coeficientes (m x m)
  ##   b    vector de términos independientes (m x 1)
  ##
  ## Retorna:
  ##   x    vector solución (1 x m)

  % Calcular triangulo superior
  [At, bt] = triang_sup(A, b);

  % Obtener solucion de la ecuacion
  x = sust_atras(At, bt);

endfunction



% ---- FUNC. AUX ----
function [At, bt] = triang_sup(A, b);
  ## Transforma un sistema de ecuaciones Ax = b en un sistema equivalente
  ## At*x = bt donde At es una matriz triangular superior
  ##
  ## Parámetros:
  ##   A    matriz de coeficientes original (m x m)
  ##   b    vector de términos independientes original (m x 1)
  ##
  ## Retorna:
  ##   At   matriz triangular superior (m x m)
  ##   bt   vector de términos independientes transformado (m x 1)

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
function x = sust_atras(A, b);
  ## Resuelve un sistema de ecuaciones Ax = b donde A es una matriz
  ## triangular superior, utilizando el método de sustitución hacia atrás
  ##
  ## Parámetros:
  ##   A    matriz triangular superior (m x m)
  ##   b    vector de términos independientes (m x 1)
  ##
  ## Retorna:
  ##   x    vector solución del sistema (1 x m)
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

