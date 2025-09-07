
clc; clear all;

% ---- FUNC. MAIN ---
% x=elimi gauss(A,b), que resuelve un sistema de ecuaciones Ax = b, utilizando las funciones
% sust atras y triang sup
% function x = elimi_gauss(A, b);



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

  bt = At(:, m + 1)
  At = At(:, 1:m)

endfunction


% ---- FUNC. AUX ----
% x=sust atras(A,b), que resuelve un sistema de ecuaciones Ax = b, donde A es triangular
% inferior, usando el metodo de sustitucion hacia atras.
% function x = sust_atras(A, b);

