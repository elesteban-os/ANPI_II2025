import numpy as np
from math import atan

def angulo(A, i, j):
    if (i != j):
        if (np.abs(A[i, j]) > 1e-16):
            return np.arctan2((2 * A[i, j]), (A[i, i] - A[j, j])) * 0.5
        else:
            return 0
    else:
        #print("Angulo: i = j")
        return None
    
def matrix_rotation(i, j, n, theta):
    Z = np.zeros((n, n))
    c = np.cos(theta)
    s = np.sin(theta)

    Z[i, i] = c - 1
    Z[j, j] = c - 1
    Z[i, j] = -s
    Z[j, i] = s

    return np.eye(n) + Z

def jacobi_valores_propios(A, iterMax, tol):

    Ak1 = A.copy()
    Ak = np.zeros_like(A)
    m = Ak.shape[0]
    ekf = lambda xk, xk1: np.linalg.norm(xk1 - xk, 2)
    xk1 = np.diag(Ak1)
    xk = np.zeros(m)
    ek = 1
    k = 0

    for k in range(1, iterMax):
        Bk = Ak1.copy()
        for i in range(0, m):
            for j in range(0, m):
                theta = angulo(Bk, i, j)
                if theta is None:
                    continue
                G = matrix_rotation(i, j, m, theta)
                Gt = G.T
                Bk = Gt @ Bk @ G

        Ak = Ak1.copy()
        Ak1 = Bk.copy()
        xk = xk1.copy()
        xk1 = np.diag(Ak1).copy()


        ek = ekf(xk, xk1)
        if (ek < tol):
            break

    return xk1, ek, k

