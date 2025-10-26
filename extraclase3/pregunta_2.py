import numpy as np
import sympy as sp
import math

A = sp.Matrix([
    [-1,  4,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0],
    [ 4, -3,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0],
    [ 0,  1,  3, -2,  0,  0,  0,  0,  0,  0,  0,  0],
    [ 0,  0, -2, -2,  3,  0,  0,  0,  0,  0,  0,  0],
    [ 0,  0,  0,  3, -1, -3,  0,  0,  0,  0,  0,  0],
    [ 0,  0,  0,  0, -3,  0,  2,  0,  0,  0,  0,  0],
    [ 0,  0,  0,  0,  0,  2,  0, -2,  0,  0,  0,  0],
    [ 0,  0,  0,  0,  0,  0, -2, -3, -2,  0,  0,  0],
    [ 0,  0,  0,  0,  0,  0,  0, -2, -2,  1,  0,  0],
    [ 0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  3,  0],
    [ 0,  0,  0,  0,  0,  0,  0,  0,  0,  3,  0, -1],
    [ 0,  0,  0,  0,  0,  0,  0,  0,  0,  0, -1,  1]
])

#----- A. Metodo de recurrencia de Sturm para polinomios caracteristicos------
# Símbolo para el parámetro .
g = sp.symbols('g')
#Recurrencia de Sturm formulas:

# p_0 = 1
# p_1 = a_11 - g
# p_k = (a_kk - g) * p_{k-1} * a_{k,k-1}^2 * p_{k-2}
def sturm(A):
    """
    Calcula los polinomios característicos p_k para k=0,...,n de la matriz simétrica A usando la recurrencia de Sturm.
    Parámetros:
    - A: matriz simétrica tridiagonal  
    Retorna:
    - polys: lista de polinomios [p_0, p_1, ..., p_n].
    """
    p0 = sp.Integer(1)
    p1 = A[0, 0] - g

    polys = [p0, p1]

    pk_2 = p0
    pk_1 = p1

    n = A.shape[0]
    for k in range(1, n):
        # expand para simplificar la expresion 
        pk = sp.expand((A[k, k] - g) * pk_1 - (A[k, k - 1] ** 2) * pk_2)
        polys.append(pk)
        pk_2, pk_1 = pk_1, pk
    return polys

polys = sturm(A)

# El polinomio característico completo (determinante) sería p_n (último elemento):
print('\nPolinomio caracteristico (p_n):')
pm = polys[-1]
print(pm)

#------B. Metodo de Gerschgorin para estimar un unico intervalo que contenga todos su valores propios ------



def gerschgorin_intervals(A):
    """
    Calcula los intervalos de Gerschgorin para una matriz simétrica A.
    Parámetros:
    - A: matriz simétrica (numpy array o sympy Matrix).
    Retorna:
    - I: lista de intervalos [min, max] para cada fila de A.
    """
    Bo = 0
    Bm = 0
    I = []
    Bj_1 = Bo
    n = A.shape[0]
    for j in range(0,n):
        if j>= n-1:
            Bj = Bm
        else:
            Bj = A[j,j+1]
        Rj = abs(Bj)+ abs(Bj_1)
        X= [A[j,j]-Rj, A[j,j]+Rj]
        I.append(X)
        Bj_1 = Bj

    return I

I = gerschgorin_intervals(A)   
# Unir todos los intervalos en uno solo que los contenga a todos
I_min = min([interval[0] for interval in I])
I_max = max([interval[1] for interval in I])
I = [I_min, I_max]
print("\nIntervalo unico que contiene todos los valores propios:", I)


#-----C. Metodo Falsa posicion para aproximar todos los valores propios -----


def subintervalos(f, a, b, n_sub=200):
    
    xs = [a + (b - a) * i / n_sub for i in range(n_sub + 1)]
    subint = []
    fa = float(f(xs[0]))
    for i in range(1, len(xs)):
        fb = float(f(xs[i]))
        if math.isnan(fa) or math.isnan(fb):
            # saltar valores no numéricos
            fa = fb
            continue
        if fa * fb <= 0:        # <=0 captura cruces y ceros
            # asegurar que no sea el mismo punto repetido
            if abs(xs[i] - xs[i-1]) > 0:
                subint.append((xs[i-1], xs[i]))
        fa = fb
    return subint


def regula_falsi(f, a, b, tol=1e-12, maxiter=200, try_subdivision=True):
  
    """
    Función que implementa el método de la falsa posición (regula falsi) para encontrar una raíz de la función f en el intervalo [a, b].
    Parámetros:
    - f: función para la cual se busca la raíz.
    - a, b: extremos del intervalo inicial.
    - tol: tolerancia para la convergencia.
    - maxiter: número máximo de iteraciones permitidas.

    Retorna:
    - c: aproximación de la raíz encontrada.
    
    """ 
    fa = float(f(a))
    fb = float(f(b))

    # caso trivial: extremo es raíz
    if abs(fa) < tol:
        return a
    if abs(fb) < tol:
        return b

    if fa * fb > 0:
        if not try_subdivision:
            raise ValueError("Los extremos no encierran una raiz y no se buscó subdivisión.")
        # buscar subintervalos con cambio de signo por muestreo
        subint = subintervalos(f, a, b, n_sub=400)
        if not subint:
            raise ValueError("No se encontró subintervalo con cambio de signo en la subdivisión.")
        # elegir el primer subintervalo
        a, b = subint[0]
        fa = float(f(a))
        fb = float(f(b))
        # si alguno extremo es raíz exacta:
        if abs(fa) < tol:
            return a
        if abs(fb) < tol:
            return b
    for _ in range(maxiter):
        c = (a * fb - b * fa) / (fb - fa)
        fc = float(f(c))
        if abs(fc) < tol or abs(b - a) < tol:
            return c
        if fa * fc < 0:
            b, fb = c, fc
        else:
            a, fa = c, fc
    return c

f_s = sp.sympify(pm)
f = sp.lambdify(g, f_s, 'numpy')
kla = regula_falsi(f, I[0], I[1])
print("\nRaíz aproximada en [", I[0], ",", I[1], "]:", kla)
