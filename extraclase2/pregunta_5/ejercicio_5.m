clc; clear; close all
Matriz = [4 -2 1; 20 -7 12; -8 13 17];

[L,U] = fact_LU(Matriz)
d = detfactlu(Matriz)
