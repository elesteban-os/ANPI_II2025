clc; clear; close all;
pkg load symbolic

syms x
resultado = unicidad_pf('((x - 3)*exp(x-3) + 2)/2', 0, 7/3);
disp(resultado);

# ---------------- FUNCIONES ---------------- #

function c = unicidad_pf(f, a, b)
## Verifica la unicidad de la solución para un sistema de ecuaciones lineales Ax = b
##
## Parámetros:
##   A    matriz de coeficientes (n x m)
##
## Retorna:
##   true   si el sistema tiene una única solución
##   false  si tiene múltiples soluciones o ninguna
% Comprueba si rango(A) es igual al número de incógnitas

  syms x
  f_sym = sym(f);
  f_d = diff(f_sym,x);

  % Paso 1: demostrar existencia
  d = existencia_pf(f_sym, f_d , a, b);

  % Paso 2: segunda derivada
  f_dd = diff(f_d, x);

  if !d
    printf('La funcion no tiene existencia \n')
    c = 0;
    return;
  else
    printf('La funcion tiene existencia \n')
  endif

  % Paso 3: encontrar puntos críticos de f'
  f_d_num = function_handle(f_d);
  puntos_criticos = solve(f_dd == 0, x);
  puntos_reales = [];

  for i = 1:length(puntos_criticos)
    punto = double(puntos_criticos(i));
    if isreal(punto) && punto >= a && punto <= b
      puntos_reales = [puntos_reales, punto];
    endif
  endfor

  % Paso 4: evaluar f' en puntos críticos y extremos
  puntos_evaluados = [f_d_num(puntos_reales), f_d_num(a), f_d_num(b)];

  % Paso 5: verificar que derivada esté en (-1,1)
  ydmax = max(puntos_evaluados);
  ydmin = min(puntos_evaluados);

  if -1 < ydmin && ydmax < 1
    printf('La funcion tiene un punto fijo \n')
    c = 1;
  else
    c = 0;
  endif
endfunction


function d = existencia_pf(f_sym, f_d, a, b)
## Verifica la existencia de solución para un sistema de ecuaciones lineales Ax = b
 ##
## Parámetros:
##   A    matriz de coeficientes (n x m)
##   b    vector de términos independientes (n x 1)
##
## Retorna:
##   true   si el sistema tiene al menos una solución
##   false  en caso contrario


  syms x
  f = function_handle(f_sym);
  try
    % Paso 2: encontrar puntos críticos de f
    puntos_criticos = solve(f_d == 0, x);
    puntos_reales = [];

    % Paso 3: filtrar reales en [a,b]
    for i = 1:length(puntos_criticos)
      punto = double(puntos_criticos(i));
      if isreal(punto) && punto >= a && punto <= b
        puntos_reales = [puntos_reales, punto];
      endif
    endfor

    % Evaluar f en puntos críticos y extremos
    puntos_evaluados = f(puntos_reales);
    if a <= f(a) && f(a) <= b
      puntos_evaluados = [puntos_evaluados, f(a)];
    endif
    if a <= f(b) && f(b) <= b
      puntos_evaluados = [puntos_evaluados, f(b)];
    endif

    ymin = min(puntos_evaluados);
    ymax = max(puntos_evaluados);

    % Paso 4: verificar intervalo
    if a <= ymin && ymax <= b
      d = 1;
    else
      d = 0;
    endif

  catch
    d = 0;
  end_try_catch
endfunction

