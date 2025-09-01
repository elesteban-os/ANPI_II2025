


function pregunta_3()
  clc; clear; close all

  % Parametros para llamar a la funcion biseccion
  f = @(x) e^x - 2 * x - 10;    % Funcion
  a = -6;                       % Intervalo
  b = -4;                       % Intervalo
  tol = 10e-8;                  % Error
  maxIter = 1000;               % Maximo de iteraciones

  % LLamada a la funcion biseccion
  [xk, k, errk] = biseccion(f, a, b, tol, maxIter)

end


function [xk k errk]=biseccion(f, a, b, tol, maxIter)

  % Obtiene la raiz de la función f en un rango [a, b]
  % por medio del metodo de la bisección
  %
  % Parametros:
  %   f        función
  %   a        limite izquierdo del rango
  %   b        limite derecho del rango
  %   tol      error máximo permitido
  %   max_iter máximo de iteraciones permitidas
  %
  % Retorna:
  %   xk     aproximación final
  %   k      cantidad de iteraciones
  %   errk   error de la función   

  % Funcion punto medio
  f_prom = @(a, b) (a + b) / 2;

  % Vectores para informacion de graficas
  err_vec = [];
  xk_vec = [];

  % Iterar la aproximación
  for k=1:maxIter
    % Calcular el punto medio entre a y b además de su valor en la funcion.
    xk = f_prom(a, b);
    fxk = f(xk);
    errk = fxk;

    % Guardar informacion para graficas
    err_vec(end + 1) = fxk;
    xk_vec(end + 1) = xk;

    % Verificar si se llego al error
    if (abs(fxk) < tol)
      break
    endif

    % Verificar teorema de Bolzano
    fa = f(a);
    if (fa * fxk) < 0
      b = xk;
    else
      a = xk;
    endif

  endfor

  % Graficas con la informacion obtenida
  figure;
  plot(1:k, abs(err_vec), 'b-o', 'LineWidth', 1.5, 'MarkerSize', 5);
  title('Iteraciones versus el error en la función');
  xlabel('Iteraciones (k)');
  ylabel('Error (err_vec)');
  grid on;

  figure;
  plot(1:k, abs(xk_vec), 'b-o', 'LineWidth', 1.5, 'MarkerSize', 5);
  title('Número de iteraciones versus la aproximación');
  xlabel('Iteraciones (k)');
  ylabel('Aproximacion (xk_vec)');
  grid on;

end



