from pregunta_1 import jacobi_valores_propios
import numpy as np

if __name__ == "__main__":
    """
    Prueba del método de Jacobi para obtener valores propios de una matriz simétrica.
    """
    m = 15
    A = np.zeros((m, m))

    for i in range(1, m + 1):
        for j in range(1, m + 1):
            A[i - 1][j - 1] = 0.5 * (i + j)

    xk, ek = jacobi_valores_propios(A, 1000, 1e-10)
    print("Valores propios aproximados:")
    print(xk)
    print("Error final:")
    print(ek)

    