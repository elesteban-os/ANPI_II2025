function d = detfactlu(A)
    # detfactlu - Calcula el determinante usando factorización LU
    # d = detfactlu(A) donde d = producto de los elementos diagonales de U

    [n, n] = size(A);


    U = fact_LU(A); #Solamente para obtener la matriz triangular superior


    d = 1;

    # Multiplicar los elementos de la diagonal
    for k = 1:n
        d = d * U(k, k);
    end
end
