clc; clear; close all
Matriz = [4 -2 1; 20 -7 12; -8 13 17];

[L,U] = fact_LU(Matriz)
d = detfactlu(Matriz)

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

