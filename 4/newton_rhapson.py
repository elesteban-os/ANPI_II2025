
from sympy import symbols, diff, lambdify, sympify
import time
f='x*exp(-x)-5-(cos(x))/(x)'
x = symbols('x')

def new_rhapson_solve(f, x0, tol=1e-8, max_iter=10000):
    """
    Solve for a root of the function f using the Newton-Raphson method.

    Parameters:"""

    if x0==0:
        raise ValueError("El valor inicial no puede ser cero.")
    
    f_sym = sympify(f, evaluate=False)
    f_num = lambdify(x, f_sym, 'numpy')
    
   
    f_d = diff(f_sym, x)
    f_d_num = lambdify(x, f_d, 'numpy')

    x_k = x0

    for i in range(max_iter):
    
        x_k1 = x_k - f_num(x_k) / f_d_num(x_k)
        #error = abs(x_k1 - x_k)
        error = abs(f_num(x_k))
        if error < tol:
            return x_k1, error, i
        x_k = x_k1
    
start = time.time()    
resultado = new_rhapson_solve(f, -0.1, tol=1e-8, max_iter=100)
end = time.time()
print("Resultado:", resultado)
print("Tiempo de ejecución:", end - start, "segundos")