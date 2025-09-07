import pandas as pd
from sympy import symbols
from sympy import symbols, lambdify, sympify, diff
import time

# Definición de variables
f='x*exp(-x)-5-(cos(x))/(x)'
x = symbols('x')
tol = 1e-8
max_iter = 10000

# Método de la bisección
def bs(f, a, b, tol=1e-8, max_iter=10000):
    """
    Obtiene la raiz de la función f en un rango [a, b]
    por medio del metodo de la bisección

    Parametros:
        f: función
        a: limite izquierdo del rango
        b: limite derecho del rango
        tol: error máximo permitido
        max_iter: máximo de iteraciones permitidas

    Retorna:
        Una tupla con los siguientes valores:
        (
            xk: aproximación final,
            fxk: error de la función,
            k: cantidad de iteraciones,
            time: tiempo de ejecución
        )  
    """

    # Funciones numericas a utilizar
    x = symbols('x')
    f_sym = sympify(f, evaluate=False)
    f_num = lambdify(x, f_sym, 'numpy')     # Función donde se aproxima el cero.
    f_prom = lambda a, b: (a + b) / 2       # Función que obtiene el punto medio entre a y b.

    # Valores iniciales
    xk, fxk = 0.0, 0.0
    k = 0

    # Iterar la aproximación
    for k in range(max_iter):
        # Calcular el punto medio y el valor en la funcion f(x)
        xk = f_prom(a, b)
        fxk = f_num(xk)

        # Verificar si se llegó al error
        if (abs(fxk) < tol):
            break

        # Verificar teorema de bolzano
        fa = f_num(a)
        if ((fa * fxk) < 0):
            b = xk
        else:
            a = xk


    # Retornar valores
    return (xk, abs(fxk), k + 1)

# Obtener resultado del método de la bisección
start_bs = time.time()
result_bs = bs(f, -0.3, -0.1, tol=tol, max_iter=max_iter)
end_bs = time.time()
time_bs = end_bs - start_bs

# Método Newton-Raphson
def nr(f, x0, tol=1e-8, max_iter=10000):
    """
    Obtiene la raiz de la función f usando el método de Newton-Raphson

    Parametros:
        f: función
        x0: valor inicial
        tol: error máximo permitido
        max_iter: máximo de iteraciones permitidas

    Retorna:
        Una tupla con los siguientes valores:
        (
            xk: aproximación final,
            fxk: error de la función,
            k: cantidad de iteraciones
        )
    
    """
    if x0==0:
        raise ValueError("El valor inicial no puede ser cero.")
    
    f_sym = sympify(f, evaluate=False)
    f_num = lambdify(x, f_sym, 'numpy')
    
   
    f_d = diff(f_sym, x)
    f_d_num = lambdify(x, f_d, 'numpy')

    x_k = x0

    for i in range(max_iter):
    
        x_k1 = x_k - f_num(x_k) / f_d_num(x_k)
        error = abs(f_num(x_k))
        if error < tol:
            return x_k1, error, i
        x_k = x_k1

# Obtener resultado del método de Newton-Raphson
start_nr = time.time()
result_nr = nr(f, -0.1, tol=tol, max_iter=max_iter)
end_nr = time.time()
time_nr = end_nr - start_nr

# Método de Steffensen
def stf(f, x0, tol=1e-8, max_iter=10000):
    """
    Función para calcular la raíz de una función con el método
    de Steffensen

    Parametros:
        f: función
        x0: valor inicial
        tol: error máximo permitido
        max_iter: máximo de iteraciones permitidas

    Retorna:
        Una tupla con los siguientes valores:
        (
            xk: aproximación final,
            fxk: error de la función,
            k: cantidad de iteraciones,
            time: tiempo de ejecución
        )  
    """

    # Funciones numericas a utilizar
    x = symbols('x')
    f_sym = sympify(f, evaluate=False)
    f_num = lambdify(x, f_sym, 'numpy')     # Función donde se aproxima el cero.
    form_steffensen = lambda xk: xk - ((f_num(xk) ** 2) / (f_num(xk + f_num(xk)) - (f_num(xk))))       # Función de la formula de Steffensen

    # Valores iniciales
    xk, fxk = x0, f_num(x0)
    k = 0

    # Restriccion
    if (f_num(xk + f_num(xk)) - f_num(xk)) != 0:
        # Iterar la aproximación
        for k in range(max_iter):
            # Calcular nuevo fxk
            fxk = f_num(xk)
            
            # Verificar si se llegó al error
            if (abs(fxk) < tol):
                break

            # Calcular nuevo xk
            xk = form_steffensen(xk)

    else:
        print("El valor inicial no es permitido")

    # Retornar valores
    return (xk, abs(fxk), k + 1)

start_stf = time.time()
result_stf = stf(f, -0.1, tol=tol, max_iter=max_iter)
end_stf = time.time()
time_stf = end_stf - start_stf

# Método Secante
def sc(f, x0, x1, tol=1e-8, max_iter=10000):
    """ 
    Función para calcular la raíz de una función con el método
    de la secante

    Parametros:
        f: función
        x0, x1: valores iniciales
        tol: error máximo permitido
        max_iter: máximo de iteraciones permitidas

    Retorna:
        Una tupla con los siguientes valores:
        (
            xk: aproximación final,
            fxk: error de la función,
            k: cantidad de iteraciones
        )    
    """

    f_sym = sympify(f, evaluate=False)
    f_num = lambdify(x, f_sym, 'numpy')
    if f_num(x1) - f_num(x0) == 0:
        x2 = "NA"
        error = "NA"
        i = "NA"
        print("La función no puede evaluarse en los puntos iniciales dados.")
        return x2, error, i

    for i in range(max_iter):
        x2 = x1 - f_num(x1)*(x1-x0)/(f_num(x1)-f_num(x0))
        error = abs(f_num(x2))
        if error < tol:
            return x2, error, i
        x0, x1 = x1, x2

# Obtener resultado del método de la secante
start_sc = time.time()
result_sc = sc(f, -0.1, -0.3, tol=tol, max_iter=max_iter)
end_sc = time.time()
time_sc = end_sc - start_sc

# Método de la Falsa Posición
def fp(f, a, b, tol=1e-8, max_iter=10000):
    """
    Obtiene la raiz de la función f en un rango [a, b]
    por medio del metodo de falsa posicion

    Parametros:
        f: función
        a: limite izquierdo del rango
        b: limite derecho del rango
        tol: error máximo permitido
        max_iter: máximo de iteraciones permitidas

    Retorna:
        Una tupla con los siguientes valores:
        (
            xk: aproximación final,
            fxk: error de la función,
            k: cantidad de iteraciones,
            time: tiempo de ejecución
        )  
    """

    # Funciones numericas a utilizar
    x = symbols('x')
    f_sym = sympify(f, evaluate=False)
    f_num = lambdify(x, f_sym, 'numpy')     # Función donde se aproxima el cero.
    f_falpos = lambda a, b: a - ((f_num(a) * (a - b)) / (f_num(a) - f_num(b)))     # Función que obtiene secante de a y b

    # Valores iniciales
    xk, fxk = 0.0, 0.0
    k = 0

    # Restriccion
    if (f_num(a) * f_num(b)) < 0:
        # Iterar la aproximación
        for k in range(max_iter):
            # Calcular el valor de xk y fxk
            xk = f_falpos(a, b)
            fxk = f_num(xk)

            # Verificar si se llegó al error
            if (abs(fxk) < tol):
                break

            # Verificar teorema de bolzano
            fa = f_num(a)
            if ((fa * fxk) < 0):
                b = xk
            else:
                a = xk
    else:
        print("Error: no cumple teorema de Bolzano")

    # Retornar valores
    return (xk, abs(fxk), k + 1)

# Obtener resultado del método de la falsa posición
start_fp = time.time()
result_fp = fp(f, -0.3, -0.1, tol=tol, max_iter=max_iter)
end_fp = time.time()
time_fp = end_fp - start_fp

# Método NHO
def hno(f, x0, tol=1e-8, max_iter=10000):
    """ 
    Método de la Hibridación de Newton y la Secante (HNO) para encontrar raíces de una función.

    Parametros:
        f : La función para la cual se busca la raíz.
        x0 : El primer valor inicial para el método.
        tol : La tolerancia para la convergencia (default es 1e-8).
        max_iter :El número máximo de iteraciones permitidas (default es 10000).
    
    Retorna:
        Una tupla con los siguientes valores:
        (
            xk: aproximación final,
            fxk: error de la función,
            k: cantidad de iteraciones
        )      
    """ 

    f_sym = sympify(f, evaluate=False)
    f_num = lambdify(x, f_sym, 'numpy')
    
   
    f_d = diff(f_sym, x)
    f_d_num = lambdify(x, f_d, 'numpy')

    xn = x0
    for i in range(max_iter):
        z_n = xn - f_num(xn) / f_d_num(xn)
        H_n = H(xn, z_n, f_d_num)
        xn1 = z_n - H_n * f_num(z_n) / f_d_num(z_n)
        error = abs(f_num(xn1))
        if error < tol:
            return xn1, error, i
        xn = xn1

# Función auxiliar para el método HNO
def H(x,z,f_d_num): 
    return (f_d_num(x)-f_d_num(z))/(3*f_d_num(z)-f_d_num(x))

# Obtener resultado del método HNO
start_hno = time.time()
result_hno = hno(f, -0.1,tol=tol, max_iter=max_iter)
end_hno = time.time()
time_hno = end_hno - start_hno

# Representacion tabulada de los resultados
data = {
    "Método": ["Bisección", "Newton-Raphson", "Steffensen", "Secante", "Falsa Posición", "HNO"],
    "Resultado Aproximado (xk)": [result_bs[0], result_nr[0], result_stf[0], result_sc[0], result_fp[0], result_hno[0]],
    "Error (|f(xk)|)": [result_bs[1], result_nr[1], result_stf[1], result_sc[1], result_fp[1], result_hno[1]],
    "Iteraciones (k)": [result_bs[2], result_nr[2], result_stf[2], result_sc[2], result_fp[2], result_hno[2]],
    "Tiempo de Ejecución (s)": [time_bs, time_nr, time_stf, time_sc, time_fp, time_hno]
}

df = pd.DataFrame(data)
print(df)

