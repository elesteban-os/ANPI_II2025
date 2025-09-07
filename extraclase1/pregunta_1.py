import numpy as np


def pentadiagonal(m, a, b, c, d, e):
    """
    Método pentadiagonal que realiza una matriz de tamaño (mxm) con los valores de 
    los vectores a, b, c, d, e.

    Parámetros:
        m: número entero positivo >= 5
        a,e: vectores tamaño mx1
        b,c: vectores tamaño (m-1)x1
        d: vector tamaño (m-2)x1

    Retorna:
        matriz: una matriz pentadiagonal
    """

    # Validación que las entradas sean correctas.
    if validar_pentadiagonal(m, a, b, c, d, e):
        matriz = np.zeros((m, m))   # Matriz donde se guarda el resultado
        
        # Llenar con vector a
        for i in range(m):
            matriz[i][i] = a[i]

        # Llenar con vector c
        for i in range(m - 1):
            matriz[i][i + 1] = c[i]
        
        # Llenar con vector e
        for i in range(m - 2):
            matriz[i][i + 2] = e[i]

        # Llenar con vector b
        for i in range(m - 1):
            matriz[i + 1][i] = b[i]
        
        # Llenar con vector d
        for i in range(m - 2):
            matriz[i + 2][i] = d[i]

        # Imprimir pentadiagonal
        return matriz
    else:
        print("Entrada invalida")

def validar_pentadiagonal(m, a, b, c, d, e):
    # Función que valida que las entradas sean correctas
    if (m >= 5 and                                   # m: un valor mayor a 5
        a.size == m and                              # a: vector tamano mx1
        b.size == m - 1 and c.size == m - 1 and      # b, c: vectores tamano (m-1)x1
        d.size == m - 2 and e.size == m - 2):        # d, e: vectores tamano (m-1)x1                       
        return True
    else:
        return False


def sistema_ecuaciones(m = 2500):
    """
    Función sistema_ecuaciones donde se construye una matriz pentadiagonal con m = 2500
    para resolver un sistema de ecuaciones.

    Imprime el resultado del error de la solución para el sistema de ecuaciones.
    """

    # Crear los vectores (se le suma 1 al rango de m)
    a = np.array([2*(i + 1) for i in range(m)])     # a(i) = 2*(i+1)
    b = np.array([(i + 1)/3 for i in range(m - 1)]) # b(i) = (i + 1)/3
    c = np.array([i/3 for i in range(m - 1)])       # c(i) = i/3
    d = np.array([(i + 2)/4 for i in range(m - 2)]) # d(i) = (i + 2)/4
    e = np.array([i/4 for i in range(m - 2)])       # e(i) = i/4

    h = np.array([2*i for i in range(m)])           # h(i) = 2*i

    # Crear la matriz pentadiagonal
    A = pentadiagonal(m, a, b, c, d, e)       

    # Solucion sistema ecuaciones
    x = np.linalg.solve(A, h)

    # Solucion ecuacion matricial
    ec_rest = np.dot(A, x) - h

    # Error con la norma euclidiana
    error = np.linalg.norm(ec_rest, 2)
    print("Sistema de ecuaciones pentadiagonal con m =", m)
    print("Error: ", error)


if __name__ == "__main__":
    sistema_ecuaciones()