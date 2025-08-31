

from sympy import symbols, lambdify, sympify
import time

def metodo_falsa_posicion(f, a, b, tol=1e-8, max_iter=10000):
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

    # Obtener el tiempo de ejecución
    time_st = time.time()

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

    # Calcular el tiempo de ejecución
    time_end = time.time() - time_st

    # Retornar valores
    return (xk, abs(fxk), k + 1, time_end)

# Función para ejecutar el método de bisección y mostrar resultados
def calc_biseccion():
    f='x*exp(-x)-5-(cos(x))/(x)'
    resultado = metodo_falsa_posicion(f, -0.3, -0.1)
    print("Aproximación:", resultado[0])
    print("Error de la función:", resultado[1])
    print("Cantidad de iteraciones:", resultado[2])
    print("Tiempo de ejecución:", resultado[3], "s")

if __name__ == "__main__":
    calc_biseccion()
