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


