import numpy as np 
import pandas as pd



def qr_gram_schmidt(A):
    """
    Realiza la descomposición QR de una matriz A utilizando el proceso de Gram-Schmidt.
    parametros:
        A    matriz de coeficientes (n x m)
    Retorna:
        Q    matriz ortogonal (n x m)
        R    matriz triangular superior (m x m)
    """
    
    A = np.array(A, dtype=float)
    n, m = A.shape

    # e_k serán las columnas de Q
    E = np.zeros((n, m))
    R = np.zeros((m, m))

    # Primer vector
    u1 = A[:, 0]
    e1 = u1 / np.linalg.norm(u1)
    E[:, 0] = e1
    R[0, 0] = np.linalg.norm(u1)

    # Resto de columnas
    for k in range(1, m):
        a_k = A[:, k]
        # u_k = a_k - sum_{j=1}^{k} <a_k, e_j> e_j
        u_k = a_k.copy()
        for j in range(k):
            R[j, k] = np.dot(a_k, E[:, j])      # <a_k, e_j>
            u_k = u_k - R[j, k] * E[:, j]       # resta proyecciones

        R[k, k] = np.linalg.norm(u_k)
        if R[k, k] == 0:
            raise ValueError("Columnas linealmente dependientes")
        e_k = u_k / R[k, k]
        E[:, k] = e_k

    Q = E
    return Q, R






A = np.array([
    [  2, -1,  3,  0,  1, -2,  4, -3,  5,  6],
    [ -3,  5, -2,  1,  4, -6,  7, -8,  9, 10],
    [  1, -2,  4, -3,  5, -7,  8, -9, 10, 11],
    [  4, -3,  1,  2, -5,  6, -7,  8, -9, 10],
    [ -5,  7, -4,  6,  3, -1,  2, -8,  9,-10],
    [  6, -4,  5, -7,  8,  1, -3,  2, -9, 10],
    [ -7,  9, -6,  8,-10, 11,  3, -5,  4,-12],
    [  8, -6,  7, -9, 10,-12, 13,  5, -4,  3],
    [ -9, 11, -8, 10,-12, 14,-16, 17,  6, -5],
    [ 10, -8,  9,-11, 12,-14, 15,-17, 18,  7]
], dtype=float)

Q, R = qr_gram_schmidt(A)
print("R=\n",pd.DataFrame(R).to_string(index=False, header=False))
print("Q=\n",pd.DataFrame(Q).to_string(index=False, header=False))