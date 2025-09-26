
function [L, U] = fact_LU(A)
  # fact_LU - Calculo de las matrices L U de A=  L U
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
   L = eye(n);  # Identidad

    for k = 1:n-1
        if abs(U(k,k)) < 1e-12
            error('Pivote numéricamente cero en (%d,%d)', k, k);
        end

        for i = k+1:n
            factor = U(i,k) / U(k,k);
            L(i,k) = factor;  # Guardar multiplicador
            U(i,k:n) = U(i,k:n) - factor * U(k,k:n);
        end
    end


endfunction





function d = subm_check(A)
  #subm_check - Verificacion de determinantes de submatrices
  n = size(A,1);
  d = 1;
  if A(1,1) == 0 # Se verificar primero A(i,j)
    d = 0;
    return;
  elseif n <= 1
    d = 1;
    return;
  end
  for k = 2:n # Se verifican las demas submatrices
    Ak = A(1:k, 1:k);
    if det(Ak) == 0
      d = 0;
      return;
    end
  end
endfunction


function g = diag_check(A,n)
  #diag_check - Verificar que la diagonal no tenga ceros (Solucion unica L U)
  for k =1:n
    if A(k,k) == 0
      g = 0;
      return;
    end
  endfor
  g = 1;
endfunction
