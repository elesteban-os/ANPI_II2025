from sympy import symbols, diff, lambdify, sympify
import time
f='x*exp(-x)-5-(cos(x))/(x)'
x = symbols('x')



def secante(f,x0, x1, tol=1e-8, max_iter=10000):
    """ 
    
    
    |"""

    f_sym = sympify(f, evaluate=False)
    f_num = lambdify(x, f_sym, 'numpy')
    if f_num(x1) - f_num(x0) == 0:
        raise ValueError("La función no puede evaluarse en los puntos iniciales dados.")
    

    for i in range(max_iter):
        x2 = x1 - f_num(x1)*(x1-x0)/(f_num(x1)-f_num(x0))
        error = abs(f_num(x2))
        if error < tol:
            return x2, error, i
        x0, x1 = x1, x2


start = time.time()
resultado = secante(f, -0.1, -0.3, tol=1e-8, max_iter=100)
end = time.time()
print("Resultado:", resultado)
print("Tiempo de ejecución:", end - start, "segundos")