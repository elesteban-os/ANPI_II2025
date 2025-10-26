from pregunta_1 import jacobi_valores_propios
import numpy as np

if __name__ == "__main__":
    # Definir A
    m = 15
    A = np.zeros((m, m))

    for i in range(1, m + 1):
        for j in range(1, m + 1):
            A[i - 1][j - 1] = 0.5 * (i + j)

    vals = jacobi_valores_propios(A, 1000, 1e-10)
    print(vals)

    eigh, _ = np.linalg.eigh(A)
    print(eigh)