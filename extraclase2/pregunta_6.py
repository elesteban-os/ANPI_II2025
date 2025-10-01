import numpy as np
import pandas as pd
def thomas(A,b): 
    """
    Resuelve un sistema de ecuaciones lineales Ax = b donde A es una matriz tridiagonal
    utilizando el algoritmo de Thomas.
    
    parametros:
        A    matriz tridiagonal de coeficientes (n x n)
        b    vector de términos independientes (n x 1)
        
    Retorna:
        x    vector solución (n x 1)
    """    

    if(check_tridiagonal(A) == False):
         print("Matriz no es tridiagonal")
         return -1

    #calcular 2 vectores pi y qi
    m ,n = A.shape
    p = np.array(np.zeros(n-1),dtype=float)
    q = np.array(np.zeros(n),dtype=float)
    
    # si i = 1, la primera iteracion
    p[0] = (A[0][1])/(A[0][0]) 
    q[0] = b[0]/A[0][0] 
    # Si i=2,...,n-1
    for i in range(1,n-1):
       k = i - 1
       p[i] = (A[i][i+1])/(A[i][i]-(p[i-1])*(A[k+1][k]))
     # Si i=2,...,n
    for i in range(1,n):
        k = i - 1
        q[i] = ((b[i]-q[i-1]*A[k+1][k])/(b[i]-p[i-1]*A[k+1][k]))
        x = np.array(np.zeros(n),dtype=float)
        x[n-1] = q[n-1]
     # Para  i = n-1,...,1
    for i in range(n-2,-1,-1):
         x[i] = q[i] - p[i]*x[i+1]

    return x



def check_tridiagonal(A):
    """
    Verifica si una matriz es tridiagonal.
    parametros:
        A    matriz de coeficientes (n x n)
    Retorna:
        True   si la matriz es tridiagonal
        False  en caso contrario
    """

    m, n = A.shape
    if m != n:
        print("Matriz no es cuadrada")
        return False
    for i in range(n): # Verificar elementos fuera de la diagonal
        for j in range(n):
            if abs(i - j) > 1 and A[i][j] != 0:
                return False
            if abs(i - j) <= 1 and A[i][j] == 0:
                return False
    return True


A= np.array([[5, 1, 0, 0], 
                   [1, 5, 1, 0],
                   [0, 1, 5, 1], 
                   [0, 0, 1, 5]], dtype=float)

b = np.array([-12,-14,-14,-12],dtype=float)

x = thomas(A,b)
print(x)
