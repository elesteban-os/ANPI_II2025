

function biseccion()

  clc; clear; close all

  % Intervalos
  a = -6; b = -4;   % [-6, -4]

  err = 10e-8;      % Error
  maxIter = 1000;   % Maximo de iteraciones
  x = (a + b) / 2;  % Valor medio entre intervalos

  # Vectores para informacion de graficas
  err_vec = [];
  xk_vec = [];

  for iter=1:maxIter
    x = (a + b) / 2;
    fx = e^x - 2 * x - 10;

    % Guardar informacion para graficas
    err_vec(end + 1) = fx;
    xk_vec(end + 1) = x;

    % Verificar si se llego al error
    if (abs(fx) < err)
      break
    endif

    % Verificar rangos
    fa = e^a - 2 * a - 10;
    if (fa * fx) < 0
      b = x;
    else
      a = x;
    endif

  endfor

  % Retornar informacion:
  [x fx iter]
  err_vec
  xk_vec


end



