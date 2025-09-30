import numpy as np
import math
import matplotlib.pyplot as plt


def pseudoinversa(matrixA, p, tol, maxIter=500):
    """
    Calcula la pseudoinversa de Moore-Penrose de una matriz utilizando
    un metodo iterativo.
    
    Parametros:
        matrixA  matriz de entrada (m x n)
        p        orden de la aproximación en la serie
        tol      tolerancia para el criterio de convergencia
        maxIter  máximo número de iteraciones permitidas
    
    Retorna:
        iters    número de iteraciones realizadas hasta convergencia
    """

    # Definir valores iniciales
    matrizAK = matrixA.T / (np.linalg.norm(matrixA, "fro") ** 2)
    m, n = matrixA.shape

    # Funcion para obtener el valor real de la sumatoria
    real_part = lambda q, p: (((-1) ** (q - 1)) * math.factorial(p)) / (math.factorial(q) * math.factorial(p - q))

    # Definir algoritmo
    iters = 0
    error = 1
    for iters in range(1, maxIter + 1):
        # Definir sumatoria
        sum_mat = np.zeros((m, m))
        for q in range(1, p + 1):
            M = matrixA @ matrizAK
            coef = real_part(q, p)
            sum_mat +=  coef * np.linalg.matrix_power(M, q - 1)
            
        # Actualizar matriz
        matrizAK = matrizAK @ sum_mat

        # Calcular error
        error = np.linalg.norm(matrixA @ matrizAK @ matrixA - matrixA, "fro")
        if error < tol:
            break
    
    # Retorn numero de iteraciones
    return iters



if __name__ == "__main__":
    # Construcción de la matriz A
    m, n = 45, 30
    A = np.zeros((m, n))
    for i in range(m):
        for j in range(n):
            A[i, j] = (i+1)**2 + (j+1)**2  

    # Diferentes valores de p
    p_vals = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]

    # Para guardar numero de iteraciones
    iteraciones = []

    # Iterar sobre diferentes valores de p
    for p in p_vals:
        iter_count = pseudoinversa(A, p, 1e-5, 100)
        iteraciones.append(iter_count)

    # Graficar resultados
    plt.plot(p_vals, iteraciones, marker='o')
    plt.xlabel('Valor de p')
    plt.ylabel('Número de iteraciones')
    plt.title('Número de iteraciones vs Valor de p')
    plt.grid(True)
    plt.show()