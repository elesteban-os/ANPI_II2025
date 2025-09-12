clc; clear; close all
pkg load symbolic

syms x
resultado = unicidad_pf('((x - 3)*exp(x-3) + 2)/2', 0, 7/3);
disp(resultado);
