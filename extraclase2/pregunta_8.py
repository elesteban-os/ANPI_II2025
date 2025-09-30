import numpy as np

def sust_atras(L, b):
    """
    Resuelve un sistema de ecuaciones lineales Lx = b donde L es una matriz
    triangular superior, utilizando el método de sustitución hacia atrás
    
    Parámetros:
        L    matriz triangular superior (n x n)
        b    vector de términos independientes (n)
    
    Retorna:
        x    vector solución del sistema (n)
    """
    x = np.zeros(np.size(b))
    for i in range(np.size(b) - 1, -1, -1):
        sum = 0
        for j in range(i + 1, np.size(b)):
            sum += L[i, j] * x[j]
        x[i] = (b[i] - sum) / L[i, i]
    return x

def gauss_seidel(A, b, x0, tol, max_iter):
    """
    Resuelve un sistema de ecuaciones lineales Ax = b utilizando el
    método iterativo de Gauss-Seidel
    
    Parámetros:
        A         matriz (n x n)
        b         vector de términos independientes (n)
        x0        vector inicial de aproximación (n)
        tol       tolerancia 
        max_iter  máximo número de iteraciones permitidas
    
    Retorna:
        errk      error final 
    """
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
    d = sust_atras(M, b)

    # Iterar
    xk = x0
    yk = 0
    errk = 0
    for iter in range(1, max_iter + 1):
        yk = -L @ xk

        # Obtener zk = M^-1 * yk
        zk = sust_atras(M, yk)

        # Obtener nuevo xk
        xk = zk + d

        # Error por medio de la norma 2
        errk = np.linalg.norm(A @ xk - b, 2)

        if (errk < tol):
            break
    
    return errk

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
    errk = gauss_seidel(A, b, x0, 1e-8, 1000)

    # Imprimir error
    print("Error:", errk)
    

