import numpy as np

# m entero positivo >= 5
# a,e vectores tamano mx1
# b,c vectores tamano (m-1)x1
# d vector tamano (m-2)x1

def pentadiagonal(m, a, b, c, d, e):
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
    # print(m, a, b, c, d, e)

    if (m >= 5 and 
        a.size == m and                              # a vector tamano mx1
        b.size == m - 1 and c.size == m - 1 and      # b, c vectores tamano (m-1)x1
        d.size == m - 2 and e.size == m - 2):        # d, e vectores tamano (m-1)x1                       
        return True
    else:
        return False


def sistema_ecuaciones(m = 2500):
    # Crear los vectores
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
    print("Error: ", error)

# # Ejemplo
# A = pentadiagonal(5, 
#     np.array([1, 2, 3, 4, 5]), 
#     np.array([6, 7, 8, 9]), 
#     np.array([10, 11, 12, 13]), 
#     np.array([14, 15, 16]), 
#     np.array([17, 18, 19])
# )
# print(A)

# h = np.array([1, 2, 3, 4, 5])
# print(h)

# # Solucion sistema ecuaciones
# x = np.linalg.solve(A, h)
# print("x = ", x)

# # Producto punto
# dot = np.dot(A, x)

# # Obtener el residuo
# resd = h - dot
# print("resd = ", resd)

# # Error con la norma euclidiana
# error = np.linalg.norm(resd, 2)
# print("error = ", error)



sistema_ecuaciones()