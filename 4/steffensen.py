

from sympy import symbols, lambdify, sympify
import time

def metodo_steffensen(f, x0, tol=1e-8, max_iter=10000):
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

    # Obtener el tiempo de ejecución
    time_st = time.time()

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

    # Calcular el tiempo de ejecución
    time_end = time.time() - time_st

    # Retornar valores
    return (xk, abs(fxk), k + 1, time_end)

# Función para ejecutar el método de bisección y mostrar resultados
def calc_steffensen():
    f='x*exp(-x) - 5 - (cos(x))/(x)'
    resultado = metodo_steffensen(f, -0.1)
    print("Aproximación:", resultado[0])
    print("Error de la función:", resultado[1])
    print("Cantidad de iteraciones:", resultado[2])
    print("Tiempo de ejecución:", resultado[3], "s")

if __name__ == "__main__":
    calc_steffensen()
