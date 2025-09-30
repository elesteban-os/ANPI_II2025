import numpy as np


def gauss_seidel(A, b, x0, tol, max_iter):
    # Obtener LDU

    # L triangular inferior estricta
    L = np.tril(A, k = -1)

    # D diagonal
    D = np.diag(np.diag(A))

    # U triangular superior estricta
    U = np.triu(A, k=1)

    # Obtener M
    M = D + U

    # Obtener d = M^-1 * b
    d = np.linalg.solve(M, b)

    # Iterar
    iter = 0
    xk = x0
    yk = 0
    errk = 0
    for iter in range(1, max_iter + 1):
        yk = -L @ xk

        # Obtener zk = M^-1 * yk
        zk = np.linalg.solve(M, yk)

        # Obtener nuevo xk
        xk = zk + d

        # Error
        errk = np.linalg.norm(A @ xk - b, 2)
        print(errk)

        if (errk < tol):
            break
    
    print("Error:", errk)

if __name__ == "__main__":
    # Definir A y b
    m = 1000
    A = np.ones((m, m))
    for i in range(0, m):
        A[i][i] = 1001
    b = np.ones(m)
    b = b.T

    # Definir x0
    x0 = np.random.rand(m)
    x0 = x0.T
    
    # Ejecutar Gauss-Seidel
    gauss_seidel(A, b, x0, 1e-8, 1000)
    

