clc; clear; close all

## ================================
## Ejercicio 5 - Factorización LU
## ================================

Matriz = [4 -2 1; 20 -7 12; -8 13 17];

[L,U] = fact_LU(Matriz)






%Parte A
function [L, U] = fact_LU(A)
## Realiza la factorización LU de una matriz A
##
## Parámetros:
##   A    matriz cuadrada (n x n)
##
## Retorna:
##   L    matriz triangular inferior (n x n)
##   U    matriz triangular superior (n x n)
##
## Nota:
##   Se verifica primero que A sea cuadrada,
##   que las submatrices principales tengan determinante distinto de cero
##   y que la diagonal no tenga ceros.
  [m, n] = size(A);
  if m ~= n
     error('La matriz debe ser cuadrada');
  end
  d = subm_check(A);
  g = diag_check(A,n);
  if [d, g] == 0
    error('Determinante de submatriz = 0 o diagonal = 0\n');
    L,U = 0;
    return;
  endif

   U = A;
   L = eye(n);  % Identidad

    for k = 1:n-1
        if abs(U(k,k)) < 1e-12
            error('Pivote numéricamente cero en (%d,%d)', k, k);
        end

        for i = k+1:n
            factor = U(i,k) / U(k,k);
            L(i,k) = factor;  % Guardar multiplicador
            U(i,k:n) = U(i,k:n) - factor * U(k,k:n);
        end
    end
endfunction


function d = subm_check(A)
## Verifica que los determinantes de las submatrices principales sean distintos de 0
##
## Parámetros:
##   A    matriz cuadrada (n x n)
##
## Retorna:
##   d    1 si todas las submatrices principales tienen determinante distinto de 0
##        0 en caso contrario

  n = size(A,1);
  d = 1;
  if A(1,1) == 0
    d = 0;
    return;
  elseif n <= 1
    d = 1;
    return;
  end
  for k = 2:n
    Ak = A(1:k, 1:k);
    if det(Ak) == 0
      d = 0;
      return;
    end
  end
endfunction



function g = diag_check(A,n)
## Verifica que la diagonal de la matriz no tenga ceros
##
## Parámetros:
##   A    matriz cuadrada (n x n)
##   n    dimensión de la matriz
##
## Retorna:
##   g    1 si no hay ceros en la diagonal
##        0 si hay al menos un cero
  for k =1:n
    if A(k,k) == 0
      g = 0;
      return;
    end
  endfor
  g = 1;
endfunction

% Parte B
function d = detfactlu(A)
## Calcula el determinante de una matriz mediante factorización LU
##
## Parámetros:
##   A    matriz cuadrada (n x n)
##
## Retorna:
##   d    determinante de la matriz A
##
## Nota:
##   Se calcula como el producto de los elementos de la diagonal de U
    [n, n] = size(A);

    U = fact_LU(A); % Solo se requiere la matriz U

    d = 1;
    for k = 1:n
        d = d * U(k, k);
    end
end


% Parte C
% Construir A
A = zeros(10, 10);
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

% Obtener determinante de A
detA = detfactlu(A)

